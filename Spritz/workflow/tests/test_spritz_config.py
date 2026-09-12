"""Tests for spritz_config.py, which is the one place config normalisation lives.

Two failures it exists to prevent, both silent.

A script cannot reach the run's config by relative path: snakemake runs from `workflow/`, where
`config/config.yaml` is the packaged default baked into the container, while the run's own config is
written into the analysis directory and passed with `--configfile`. Opening the relative path parses
fine and yields the defaults, so a bacterial run reported itself as Homo_sapiens.

And `--config vcf=x.vcf`, like a hand-written `vcf: "x.vcf"`, gives a string. `check` is a length
test that any non-empty string passes, so the VCF branch was taken and then the string iterated one
character at a time.

This replaces an earlier test file that reimplemented the coercion locally and so only exercised its
own copy: reinstating the original bug in the real implementation left all of its cases passing.
"""
import os
import sys
import textwrap

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

# noqa: E402 below - the sys.path insert above has to run first
from spritz_config import check, choose_configfile, load, normalise  # noqa: E402


@pytest.mark.parametrize(
    "given,expected",
    [
        ("a.vcf", ["a.vcf"]),                      # --config vcf=a.vcf
        ("", []),                                   # the empty default written as a scalar
        (["a.vcf"], ["a.vcf"]),
        (["a.vcf", "b.vcf"], ["a.vcf", "b.vcf"]),
        ([], []),
    ],
)
def test_a_scalar_vcf_becomes_a_one_element_list(given, expected):
    assert normalise({"vcf": given})["vcf"] == expected


def test_a_coerced_scalar_yields_one_entry_not_one_per_character():
    """The actual bug: iterating "a.vcf" gives 'a', '.', 'v', 'c', 'f'."""
    assert list(normalise({"vcf": "sample_one.vcf"})["vcf"]) == ["sample_one.vcf"]


def test_an_empty_scalar_reads_as_absent_rather_than_present():
    """It must not become [""], which check() would pass for a filename that is not there."""
    config = normalise({"vcf": ""})
    assert config["vcf"] == []
    assert check(config, "vcf") is False


def test_bacteria_are_coerced_to_bootstrap():
    """Ensembl Bacteria publishes no variation, so the downloaded route can never apply."""
    assert normalise({"division": "bacteria"})["known_sites"] == "bootstrap"
    assert normalise({"division": "bacteria", "known_sites": "ensembl"})["known_sites"] == "bootstrap"


def test_the_bacteria_coercion_folds_case():
    """-e=Bacteria is accepted by the CLI, so the workflow must not read it as a vertebrate."""
    assert normalise({"division": "Bacteria"})["known_sites"] == "bootstrap"


def test_an_explicit_bootstrap_is_left_alone_for_vertebrates():
    assert normalise({"division": "vertebrates", "known_sites": "bootstrap"})["known_sites"] == "bootstrap"


def test_a_vertebrate_is_not_coerced():
    """The coercion must not quietly move a vertebrate off the downloaded route."""
    assert normalise({"division": "vertebrates", "known_sites": "ensembl"})["known_sites"] == "ensembl"
    assert "known_sites" not in normalise({"division": "vertebrates"})


def test_load_reads_the_path_it_is_given_not_the_relative_default(tmp_path, monkeypatch):
    """The bug this module exists for: the relative path is the packaged default, not the run."""
    packaged = tmp_path / "config"
    packaged.mkdir()
    (packaged / "config.yaml").write_text('species: "Homo_sapiens"\ndivision: "vertebrates"\n')
    real = tmp_path / "run.yaml"
    real.write_text('species: "Pseudomonas_aeruginosa_pao1_gca_000006765"\ndivision: "bacteria"\n')

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("SPRITZ_CONFIG", raising=False)
    assert load()["species"] == "Homo_sapiens"                  # what the old code always got
    assert load(str(real))["division"] == "bacteria"            # explicit path
    monkeypatch.setenv("SPRITZ_CONFIG", str(real))
    assert load()["division"] == "bacteria"                     # what the rules now set


def test_load_applies_the_normalisations(tmp_path, monkeypatch):
    path = tmp_path / "run.yaml"
    path.write_text(textwrap.dedent('''\
        division: "bacteria"
        vcf: "one.vcf"
        '''))
    config = load(str(path))
    assert config["vcf"] == ["one.vcf"]
    assert config["known_sites"] == "bootstrap"


def test_an_empty_config_file_does_not_crash(tmp_path):
    path = tmp_path / "empty.yaml"
    path.write_text("")
    assert load(str(path)) == {}


class TestChooseConfigfile:
    """Which entry of snakemake's workflow.configfiles is the run's own config.

    The list is seeded with the --configfile paths and then appended to by the `configfile:`
    directive at the top of the Snakefile, so its last entry is the packaged default. Taking
    [-1] shipped a fix that did nothing: a bacterial cluster run downloaded the human UniProt
    proteome under the bacterium's filename and wrote a prose.txt describing a quant run.
    """

    def test_picks_the_supplied_config_not_the_directive_default(self, tmp_path):
        default = tmp_path / "workflow" / "config" / "config.yaml"
        default.parent.mkdir(parents=True)
        default.write_text("species: Homo_sapiens\n")
        run = tmp_path / "results" / "config" / "config.yaml"
        run.parent.mkdir(parents=True)
        run.write_text("species: Pseudomonas\n")

        # The order snakemake actually produces: CLI first, directive appended last.
        assert choose_configfile([str(run), str(default)], str(default)) == str(run)

    def test_falls_back_to_the_default_when_nothing_was_supplied(self, tmp_path):
        default = tmp_path / "config" / "config.yaml"
        default.parent.mkdir(parents=True)
        default.write_text("species: Homo_sapiens\n")
        assert choose_configfile([str(default)], str(default)) == str(default)

    def test_empty_list_falls_back_to_the_default(self):
        assert choose_configfile([], "config/config.yaml") == "config/config.yaml"

    def test_recognises_the_default_through_a_different_spelling_of_the_same_path(self, tmp_path):
        """The directive's path and the one built from workflow.basedir need not match textually."""
        default = tmp_path / "config" / "config.yaml"
        default.parent.mkdir(parents=True)
        default.write_text("species: Homo_sapiens\n")
        indirect = tmp_path / "config" / ".." / "config" / "config.yaml"
        assert choose_configfile([str(indirect)], str(default)) == str(default)

    def test_last_supplied_config_wins(self, tmp_path):
        """--configfile a b means b overrides a, so b is the run's effective config."""
        default = tmp_path / "config" / "config.yaml"
        default.parent.mkdir(parents=True)
        default.write_text("species: Homo_sapiens\n")
        first, second = tmp_path / "a.yaml", tmp_path / "b.yaml"
        first.write_text("species: A\n")
        second.write_text("species: B\n")
        chosen = choose_configfile([str(first), str(second), str(default)], str(default))
        assert chosen == str(second)
