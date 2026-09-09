"""Tests the scalar-to-list coercion common.smk applies to the `vcf` config key.

The failure this prevents is quiet and total. `--config vcf=a.vcf` on the command line, and a
hand-written `vcf: "a.vcf"` in config.yaml, both give a plain string. `check()` is
`len(config[field]) > 0`, which is true for any non-empty string, so the VCF branch is taken - and
then `[posixpath.join(dir, v) for v in config['vcf']]` iterates the string's *characters* and asks
snakemake for one input file per letter.

common.smk is a snakefile and cannot be imported, so the coercion is reproduced here and pinned
against the same cases. test_smoke_commands.py sets the precedent for checking a workflow file's
content from a test rather than executing it.
"""
import os
import re

import pytest

COMMON_SMK = os.path.join(os.path.dirname(__file__), "..", "rules", "common.smk")


def coerce(value):
    """The coercion under test, mirroring common.smk."""
    config = {"vcf": value}
    if isinstance(config.get("vcf"), str):
        config["vcf"] = [config["vcf"]] if config["vcf"] else []
    return config["vcf"]


def check(value):
    """check() from common.smk."""
    return value is not None and len(value) > 0


@pytest.mark.parametrize(
    "given,expected",
    [
        ("a.vcf", ["a.vcf"]),           # --config vcf=a.vcf
        ("", []),                        # the empty default, written as a scalar
        (["a.vcf"], ["a.vcf"]),         # already a list
        (["a.vcf", "b.vcf"], ["a.vcf", "b.vcf"]),
        ([], []),
    ],
)
def test_a_scalar_vcf_becomes_a_one_element_list(given, expected):
    assert coerce(given) == expected


def test_a_coerced_scalar_yields_one_path_not_one_per_character():
    """The actual bug: iterating "a.vcf" gives 'a', '.', 'v', 'c', 'f'."""
    assert list(coerce("a.vcf")) == ["a.vcf"]
    assert len(list(coerce("sample_one.vcf"))) == 1


def test_an_empty_scalar_reads_as_absent_rather_than_present():
    """len("") is 0 so check() already said absent, but the coercion must not turn it into [""]."""
    assert check(coerce("")) is False
    assert coerce("") == []


def test_common_smk_still_performs_the_coercion():
    """A guard against the coercion being dropped, which no other test here would notice."""
    with open(COMMON_SMK, encoding="utf-8") as handle:
        source = handle.read()
    assert re.search(r'isinstance\(\s*config\.get\("vcf"\)\s*,\s*str\s*\)', source), \
        "common.smk no longer coerces a scalar vcf config value to a list"
    # And it must happen before check() is used to decide the branch.
    assert source.index('isinstance(config.get("vcf")') < source.index("def check(field)")
