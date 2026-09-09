"""Tests for ensembl_bacteria.py, which resolves an Ensembl Bacteria reference.

Three things here are not guessable and would each produce a 404 rather than a wrong answer, so they
are pinned against real values taken from Ensembl Genomes release 63:

- the collection directory, which is not derivable from the species name (the 57 Pseudomonas
  aeruginosa strains are spread over about 30 collections) and comes from column 14 of
  species_EnsemblBacteria.txt;
- the DNA filename's trailing underscore on the assembly token, which is why filenames are read from
  the per-directory CHECKSUMS index rather than templated;
- the assembly-name sanitisation Ensembl applies when it builds a filename.

Unlike the other script tests these call the functions directly rather than running the script as a
subprocess: the parts worth testing are pure, and exercising the rest would mean reaching the network.
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts")))

from ensembl_bacteria import (  # noqa: E402  - needs the path set above
    parse_metadata,
    pick,
    resolve,
    sanitise_assembly,
)

HEADER = (
    "#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\t"
    "microarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\t"
    "species_id\n"
)


def row(name, species, assembly, collection, species_id="166"):
    return (
        f"{name}\t{species}\tEnsemblBacteria\t208964\t{assembly}\tGCA_000006765.1\t2022-12-Prokka\t"
        f"N\tN\tN\tN\tN\tN\t{collection}_core_63_116_1\t{species_id}\n"
    )


# Verbatim from release 63.
PAO1 = row(
    "Pseudomonas aeruginosa PAO1 (GCA_000006765)",
    "pseudomonas_aeruginosa_pao1_gca_000006765",
    "ASM676v1",
    "bacteria_5_collection",
)
BL04 = row(
    "Pseudomonas aeruginosa BL04 (GCA_000481065)",
    "pseudomonas_aeruginosa_bl04_gca_000481065",
    "Pseu_aeru_BL04_V1",
    "bacteria_119_collection",
)

METADATA = HEADER + PAO1 + BL04


def test_the_collection_comes_from_the_core_db_column():
    """It is not derivable from the name: these two strains are in different collections."""
    parsed = parse_metadata(METADATA)
    assert parsed["pseudomonas_aeruginosa_pao1_gca_000006765"] == ("ASM676v1", "bacteria_5_collection")
    assert parsed["pseudomonas_aeruginosa_bl04_gca_000481065"] == (
        "Pseu_aeru_BL04_V1", "bacteria_119_collection")


def test_the_comment_header_is_not_read_as_a_species():
    assert "#name" not in parse_metadata(METADATA)
    assert len(parse_metadata(METADATA)) == 2


def test_a_bare_species_name_fails_with_the_strains_that_do_exist():
    """The 2020-era name. Every species directory now carries a _gca_ suffix, so it cannot match."""
    with pytest.raises(SystemExit) as caught:
        resolve(parse_metadata(METADATA), "pseudomonas_aeruginosa")
    message = str(caught.value)
    assert "no species directory named 'pseudomonas_aeruginosa'" in message
    assert "_gca_" in message
    # Names real alternatives rather than only reporting the failure.
    assert "pseudomonas_aeruginosa_pao1_gca_000006765" in message


def test_a_species_matched_nowhere_says_so_rather_than_suggesting_nothing():
    with pytest.raises(SystemExit) as caught:
        resolve(parse_metadata(METADATA), "escherichia_coli_gca_000005845")
    assert "No species directory contains that string" in str(caught.value)


@pytest.mark.parametrize(
    "assembly,expected",
    [
        # Every case here is a real assembly name from release 63.
        ("ASM676v1", "ASM676v1"),
        ("IMG-taxon 2667527408 annotated assembly", "IMG-taxon_2667527408_annotated_assembly"),
        ("17870_2#15", "17870_2_15"),
        ("BRSU_AN4859/03", "BRSU_AN4859_03"),
        ("EC_O111:H8_CVM9634_1.0", "EC_O111_H8_CVM9634_1.0"),
        ("NGEN (DNA Star)", "NGEN_DNA_Star_"),
        ("SBR5(T)", "SBR5_T_"),
    ],
)
def test_assembly_names_are_sanitised_the_way_ensembl_names_its_files(assembly, expected):
    """2,217 of the 31,332 assemblies contain a space and 1,580 contain a '#'."""
    assert sanitise_assembly(assembly) == expected


def test_dots_and_hyphens_survive_sanitisation():
    """So a filename cannot be split on '.' to recover the fields - 739 assemblies contain a dot."""
    assert sanitise_assembly("D15-8W.seq.sqn for a Marinobacter nanhaiticus assembly, version 1.0") == (
        "D15-8W.seq.sqn_for_a_Marinobacter_nanhaiticus_assembly_version_1.0")


# The real CHECKSUMS listing for PAO1's dna/ directory. Note the trailing underscore on the assembly
# token, which is present for most strains and absent for some, with nothing in the metadata to say
# which - the reason filenames are read from here rather than constructed.
PAO1_DNA_FILES = [
    "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1_.dna.nonchromosomal.fa.gz",
    "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1_.dna.toplevel.fa.gz",
    "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1_.dna_rm.toplevel.fa.gz",
    "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1_.dna_sm.toplevel.fa.gz",
]


def test_the_soft_and_repeat_masked_genomes_are_not_picked():
    """dna_rm and dna_sm end in '_rm.toplevel'/'_sm.toplevel', so '.dna.toplevel' excludes them."""
    assert pick(PAO1_DNA_FILES, ".dna.toplevel.fa.gz", "dna/") == (
        "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1_.dna.toplevel.fa.gz")


def test_a_chromosome_level_gff3_is_not_picked_over_the_whole_one():
    """E. coli K-12 ships both; only the one numbered with the release is wanted."""
    names = [
        "Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.63.chromosome.Chromosome.gff3.gz",
        "Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.63.gff3.gz",
    ]
    assert pick(names, ".63.gff3.gz", "gff3/") == (
        "Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.63.gff3.gz")


def test_the_abinitio_proteome_is_not_picked_over_the_real_one():
    names = [
        "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.pep.abinitio.fa.gz",
        "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.pep.all.fa.gz",
    ]
    assert pick(names, ".pep.all.fa.gz", "pep/") == (
        "Pseudomonas_aeruginosa_pao1_gca_000006765.ASM676v1.pep.all.fa.gz")


def test_an_ambiguous_or_missing_match_errors_rather_than_guessing():
    """Silently taking the first would download the wrong genome and never say so."""
    with pytest.raises(SystemExit) as missing:
        pick(PAO1_DNA_FILES, ".pep.all.fa.gz", "dna/")
    assert "found 0" in str(missing.value)

    with pytest.raises(SystemExit) as ambiguous:
        pick(["a.dna.toplevel.fa.gz", "b.dna.toplevel.fa.gz"], ".dna.toplevel.fa.gz", "dna/")
    assert "found 2" in str(ambiguous.value)


def test_the_taxonomy_id_is_read_from_the_right_column():
    """get_proteome.py looks a bacterial proteome up by this: the name UniProt files PAO1 under is
    "Pseudomonas aeruginosa (strain ATCC 15692 / ... / PAO1)", which no name-based search reaches."""
    from ensembl_bacteria import parse_taxonomy_ids
    ids = parse_taxonomy_ids(METADATA)
    assert ids["pseudomonas_aeruginosa_pao1_gca_000006765"] == "208964"
    assert "#name" not in ids


def test_the_taxonomy_column_is_not_confused_with_the_assembly_column():
    """They are adjacent - taxonomy_id is column 4 and assembly column 5 - so an off-by-one here
    would send a plausible-looking but wrong value to UniProt."""
    from ensembl_bacteria import parse_metadata, parse_taxonomy_ids
    species = "pseudomonas_aeruginosa_pao1_gca_000006765"
    assert parse_taxonomy_ids(METADATA)[species] == "208964"
    assert parse_metadata(METADATA)[species][0] == "ASM676v1"
