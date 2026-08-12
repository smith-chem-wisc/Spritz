"""Tests for the genomes.csv generator.

These pin the two transformations that a regeneration gets wrong by default, plus the
four-column contract that Spritz's reader depends on. No network access.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts"))

import update_genomes as ug


class TestAssemblyPatchSuffix:
    """Ensembl reports GRCh38.p14; the files it serves are named GRCh38."""

    def test_human_patch_suffix_is_removed(self):
        assert ug.strip_patch("GRCh38.p14") == "GRCh38"

    @pytest.mark.parametrize("assembly", ["GRCm39", "GRCz11", "ROS_Cfam_1.0", "Mmul_10", "BDGP6.46"])
    def test_unpatched_assemblies_are_untouched(self, assembly):
        assert ug.strip_patch(assembly) == assembly

    def test_a_trailing_version_that_is_not_a_patch_survives(self):
        # BDGP6.46 ends in .46 but is not a patch; only .p<N> is a patch.
        assert ug.strip_patch("BDGP6.46") == "BDGP6.46"


class TestCommaSanitising:
    """The reader splits on commas and requires exactly four fields."""

    def test_a_comma_in_a_species_name_is_replaced(self):
        assert ug.sanitise("Caenorhabditis elegans (Nematode, N2)") == "Caenorhabditis elegans (Nematode; N2)"

    def test_a_sanitised_row_still_has_four_fields(self):
        row = ("release-116", "caenorhabditis_elegans", ug.sanitise("C. elegans (Nematode, N2)").lower(), "WBcel235")
        assert len(",".join(row).split(",")) == 4


class TestMerge:
    def test_historical_rows_are_kept(self):
        existing = [("release-97", "homo_sapiens", "human", "GRCh38")]
        fetched = [("release-116", "homo_sapiens", "human", "GRCh38")]
        assert len(ug.build(existing, fetched)) == 2

    def test_fetched_rows_win_over_stale_ones(self):
        existing = [("release-116", "mus_musculus", "mouse", "GRCm38")]
        fetched = [("release-116", "mus_musculus", "mouse", "GRCm39")]
        assert ug.build(existing, fetched) == [("release-116", "mus_musculus", "mouse", "GRCm39")]

    def test_rows_are_ordered_numerically_not_lexically(self):
        # "release-97" sorts after "release-116" as a string, which would put the newest
        # release in the middle of the file and change what the GUI lists first.
        rows = ug.build([], [("release-97", "a", "a", "A"), ("release-116", "a", "a", "A")])
        assert [r[0] for r in rows] == ["release-97", "release-116"]


class TestGeneratedFile:
    """Checks against the committed file, so a bad regeneration cannot land unnoticed."""

    @staticmethod
    def rows():
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "SpritzBackend", "genomes.csv")
        with open(path, encoding="utf-8") as handle:
            return [line.rstrip("\n") for line in handle if line.strip() and not line.startswith("#")]

    def test_every_row_has_exactly_four_fields(self):
        bad = [r for r in self.rows() if len(r.split(",")) != 4]
        assert bad == [], f"rows Spritz would reject: {bad[:5]}"

    def test_no_assembly_carries_a_patch_suffix(self):
        bad = [r for r in self.rows() if ug.PATCH_SUFFIX.search(r.split(",")[3])]
        assert bad == [], f"these would build 404 URLs: {bad[:5]}"

    def test_the_current_human_and_mouse_assemblies_are_present(self):
        rows = set(self.rows())
        assert any(r.startswith("release-116,homo_sapiens,") and r.endswith(",GRCh38") for r in rows)
        assert any(r.startswith("release-116,mus_musculus,") and r.endswith(",GRCm39") for r in rows)

    def test_species_with_their_own_gene_set_version_are_listed(self):
        # These are numbered 63 in release 116 because Ensembl imports them from WormBase, FlyBase
        # and SGD. They are usable because the download rule discovers the real gff3 name, so
        # excluding them would drop three of the most commonly used model organisms.
        rows = self.rows()
        for species in ("caenorhabditis_elegans", "drosophila_melanogaster", "saccharomyces_cerevisiae"):
            assert any(r.startswith(f"release-116,{species},") for r in rows), species

    def test_no_release_name_is_a_prefix_of_another(self):
        # Spritz groups rows by release with a substring match, so if one release name were
        # a prefix of another (release-11 and release-116) species would be attributed to
        # the wrong release and the dictionary build would throw on the duplicate.
        releases = {r.split(",")[0] for r in self.rows()}
        clashes = [(a, b) for a in releases for b in releases if a != b and b.startswith(a)]
        assert clashes == [], f"prefix clashes: {clashes}"
