"""Tests for snpeff_config.py, which decides the codon tables SnpEff uses per chromosome.

The stakes are that SnpEff falls back to the standard table for any chromosome the config does not
name, so a missing or misspelled contig name silently mistranslates mitochondrial genes rather than
failing. These tests pin the contig name and table for each species Ensembl carries that is not a
vertebrate, and pin the case-folding, since config.yaml ships the species capitalised.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from snpeff_config import (  # noqa: E402  - needs the path set above
    KNOWN_CODON_TABLES,
    MITOCHONDRIA_BY_SPECIES,
    VERTEBRATE_MITOCHONDRIA,
    mitochondria,
    snpeff_config_block,
)


def codon_table_lines(block):
    """The codonTable lines of a block, stripped of the leading tab."""
    return [line.strip() for line in block.splitlines() if ".codonTable" in line]


def test_a_vertebrate_gets_table_2_on_mt():
    assert mitochondria("homo_sapiens") == VERTEBRATE_MITOCHONDRIA
    lines = codon_table_lines(
        snpeff_config_block("homo_sapiens.GRCh38", "homo_sapiens", "GRCh38", "111")
    )
    assert lines == [
        "homo_sapiens.GRCh38.MT.codonTable : Vertebrate_Mitochondrial",
        "homo_sapiens.GRCh38.M.codonTable : Vertebrate_Mitochondrial",
    ]


def test_an_unlisted_species_falls_back_to_the_vertebrate_default():
    """Ensembl's main site is vertebrates, so the default is right for everything not listed."""
    assert mitochondria("danio_rerio") == VERTEBRATE_MITOCHONDRIA
    assert mitochondria("gallus_gallus") == VERTEBRATE_MITOCHONDRIA


# The contig names come from the Ensembl REST assembly info for each species; the tables are the
# NCBI translation tables for the clade. Both halves matter: a wrong contig name means no line
# applies and the genes stay on Standard.
@pytest.mark.parametrize(
    "species, expected",
    [
        ("caenorhabditis_elegans", [("MtDNA", "Invertebrate_Mitochondrial")]),
        ("drosophila_melanogaster", [("mitochondrion_genome", "Invertebrate_Mitochondrial")]),
        ("saccharomyces_cerevisiae", [("Mito", "Yeast_Mitochondrial")]),
        ("ciona_intestinalis", [("MT", "Ascidian_Mitochondrial")]),
        ("ciona_savignyi", []),
    ],
)
def test_non_vertebrates_get_their_own_contig_and_table(species, expected):
    assert list(mitochondria(species)) == expected


def test_ciona_intestinalis_is_not_left_on_the_vertebrate_table():
    """This is the species the old hardcoded MT line reached and mislabelled.

    Ascidian mitochondria read AGA and AGG as glycine where vertebrates read them as stop, so the
    difference is not cosmetic.
    """
    lines = codon_table_lines(
        snpeff_config_block("ciona_intestinalis.KH", "ciona_intestinalis", "KH", "111")
    )
    assert lines == ["ciona_intestinalis.KH.MT.codonTable : Ascidian_Mitochondrial"]
    assert "Vertebrate_Mitochondrial" not in "".join(lines)


def test_worm_fly_and_yeast_are_not_left_on_the_standard_table():
    """These three name their contig something other than MT, so a hardcoded MT line missed them.

    Missing entirely is worse than wrong: the fallback is Standard, where TGA is a stop instead of
    tryptophan.
    """
    for species, assembly, contig in [
        ("caenorhabditis_elegans", "WBcel235", "MtDNA"),
        ("drosophila_melanogaster", "BDGP6.54", "mitochondrion_genome"),
        ("saccharomyces_cerevisiae", "R64-1-1", "Mito"),
    ]:
        genome = f"{species}.{assembly}"
        lines = codon_table_lines(snpeff_config_block(genome, species, assembly, "111"))
        assert len(lines) == 1
        assert lines[0].startswith(f"{genome}.{contig}.codonTable : ")


def test_species_lookup_folds_case():
    """config.yaml ships species capitalised, e.g. "Homo_sapiens"."""
    assert mitochondria("Saccharomyces_cerevisiae") == mitochondria("saccharomyces_cerevisiae")
    assert mitochondria("CAENORHABDITIS_ELEGANS") == mitochondria("caenorhabditis_elegans")
    lines = codon_table_lines(
        snpeff_config_block(
            "Saccharomyces_cerevisiae.R64-1-1", "Saccharomyces_cerevisiae", "R64-1-1", "111"
        )
    )
    assert lines == ["Saccharomyces_cerevisiae.R64-1-1.Mito.codonTable : Yeast_Mitochondrial"]


def test_every_table_emitted_is_one_snpeff_defines():
    """A table SnpEff has no codon.* entry for would fall back rather than fail."""
    tables = {table for _, table in VERTEBRATE_MITOCHONDRIA}
    for pairs in MITOCHONDRIA_BY_SPECIES.values():
        tables.update(table for _, table in pairs)
    assert tables <= KNOWN_CODON_TABLES


def test_the_block_is_shaped_the_way_snpeff_reads_it():
    block = snpeff_config_block("homo_sapiens.GRCh38", "homo_sapiens", "GRCh38", "111")
    lines = block.split("\n")
    # A blank first line separates this entry from whatever precedes it in snpEff.config.
    assert lines[0] == ""
    assert lines[1] == "# homo_sapiens.GRCh38"
    assert lines[2] == "homo_sapiens.GRCh38.genome : homo_sapiens GRCh38"
    assert lines[3] == "homo_sapiens.GRCh38.reference : https://ftp.ensembl.org/pub/release-111/"
    # codonTable entries are indented under the genome, and the block ends with a newline.
    assert all(line.startswith("\t") for line in lines[4:] if line)
    assert block.endswith("\n")


def test_the_block_contains_no_literal_backslash_n():
    """`echo "\\n# ref"` in bash wrote a literal backslash-n into the config; printf does not."""
    block = snpeff_config_block("homo_sapiens.GRCh38", "homo_sapiens", "GRCh38", "111")
    assert "\\n" not in block


def test_a_custom_gene_model_says_so_rather_than_claiming_to_be_the_reference():
    block = snpeff_config_block(
        "combined.transcripts.genome.gff3",
        "homo_sapiens",
        "GRCh38",
        "111",
        gene_model="StringTie-assembled transcripts",
    )
    assert (
        "combined.transcripts.genome.gff3.genome : homo_sapiens GRCh38 with "
        "StringTie-assembled transcripts" in block
    )
    # The custom gene model still needs the species' codon tables.
    assert codon_table_lines(block) == [
        "combined.transcripts.genome.gff3.MT.codonTable : Vertebrate_Mitochondrial",
        "combined.transcripts.genome.gff3.M.codonTable : Vertebrate_Mitochondrial",
    ]


def test_no_entry_describes_the_genome_as_human_refseq():
    """Every species used to be labelled "Human genome ... using RefSeq transcripts"."""
    block = snpeff_config_block("danio_rerio.GRCz11", "danio_rerio", "GRCz11", "111")
    assert "Human" not in block
    assert "RefSeq" not in block
    assert "ncbi.nlm.nih.gov" not in block
