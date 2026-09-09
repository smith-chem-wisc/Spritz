"""Tests for stage_user_vcf.py, which brings a user-supplied VCF into the workflow.

The stakes are a silent empty database. SnpEff calls a variant on an unknown contig
ERROR_CHROMOSOME_NOT_FOUND and exits 0, so a VCF with UCSC contig names run against an Ensembl
reference produces a protein database with no variants in it and no failure anywhere. These tests pin
that a total mismatch stops the run, and - as importantly - that a partial one does not, since a VCF
naming scaffolds the primary assembly omits is ordinary.
"""
import gzip
import textwrap

import pytest

from conftest import run_script

ENSEMBL_FAI = "1\t248956422\t112\t60\t61\n2\t242193529\t253105774\t60\t61\nMT\t16569\t0\t60\t61\n"

VCF_HEADER = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##source=someothercaller
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A
    """
)

# Sample columns are not decoration here: the database builder keeps only variants carrying
# genotypes, and reads AD per allele index, so both are preconditions rather than nice-to-haves.
CALLED = "GT:AD\t0/1:10,8"


def vcf(*data_lines):
    """Data lines given as the first eight columns; a called genotype with depths is appended."""
    return VCF_HEADER + "".join(f"{line}\t{CALLED}\n" for line in data_lines)


@pytest.fixture
def fai(tmp_path):
    path = tmp_path / "genome.fa.fai"
    path.write_text(ENSEMBL_FAI)
    return path


def stage(tmp_path, fai, text, gzipped=False, name="user.vcf"):
    path = tmp_path / (name + (".gz" if gzipped else ""))
    if gzipped:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return run_script("stage_user_vcf.py", tmp_path, args=(str(path), str(fai)))


def test_a_matching_vcf_passes_through_byte_for_byte(tmp_path, fai):
    text = vcf("1\t1000\t.\tA\tG\t50\tPASS\t.", "MT\t20\t.\tC\tT\t50\tPASS\t.")
    result = stage(tmp_path, fai, text)
    assert result.returncode == 0, result.stderr
    assert result.stdout == text
    assert "Staged 2 variant(s) on 2 reference contig(s) for 1 sample(s): SAMPLE_A." in result.stderr


def test_a_gzipped_vcf_is_decompressed(tmp_path, fai):
    # Detected by magic number, not by extension: a plain .vcf is sometimes gzipped anyway, and the
    # rest of the workflow hands this file to SnpEff and to cp as text either way.
    text = vcf("1\t1000\t.\tA\tG\t50\tPASS\t.")
    result = stage(tmp_path, fai, text, gzipped=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout == text


def test_ucsc_contig_names_stop_the_run(tmp_path, fai):
    """The failure this script exists for: every variant would be dropped, silently."""
    text = vcf("chr1\t1000\t.\tA\tG\t50\tPASS\t.", "chrM\t20\t.\tC\tT\t50\tPASS\t.")
    result = stage(tmp_path, fai, text)
    assert result.returncode != 0
    assert "no variant" in result.stderr
    # Names the way out, not just the problem.
    assert "convert_ucsc2ensembl.py" in result.stderr
    # Shows both vocabularies, so the mismatch is visible rather than merely asserted.
    assert "chr1" in result.stderr and "MT" in result.stderr


def test_a_partly_matching_vcf_is_staged_with_a_warning(tmp_path, fai):
    """A scaffold the primary assembly omits is ordinary and must not stop the run."""
    text = vcf(
        "1\t1000\t.\tA\tG\t50\tPASS\t.",
        "KI270728.1\t50\t.\tA\tT\t50\tPASS\t.",
        "KI270728.1\t90\t.\tG\tC\t50\tPASS\t.",
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode == 0, result.stderr
    assert result.stdout == text
    assert "2 variant(s) on 1 contig(s) absent from the reference" in result.stderr
    assert "KI270728.1" in result.stderr


def test_a_header_only_vcf_stops_the_run(tmp_path, fai):
    """No variants at all is a mistake worth catching, not an empty success."""
    result = stage(tmp_path, fai, VCF_HEADER)
    assert result.returncode != 0
    assert "no variant" in result.stderr


def test_an_empty_reference_index_stops_the_run(tmp_path):
    """Rather than reporting every contig as unmatched and blaming the user's VCF."""
    empty = tmp_path / "empty.fa.fai"
    empty.write_text("")
    result = stage(tmp_path, empty, vcf("1\t1000\t.\tA\tG\t50\tPASS\t."))
    assert result.returncode != 0
    assert "the .fai is empty" in result.stderr


def test_the_unmatched_contig_list_is_capped(tmp_path, fai):
    """A VCF against the wrong assembly has thousands of unknown contigs; do not print them all."""
    text = vcf(*[f"scaffold_{i}\t10\t.\tA\tG\t50\tPASS\t." for i in range(50)],
               "1\t1000\t.\tA\tG\t50\tPASS\t.")
    result = stage(tmp_path, fai, text)
    assert result.returncode == 0, result.stderr
    assert "50 variant(s) on 50 contig(s) absent" in result.stderr
    assert result.stderr.count("scaffold_") == 10


def test_a_sites_only_vcf_stops_the_run(tmp_path, fai):
    """No sample columns means no genotypes, and the builder discards variants without them - so
    this would otherwise produce a database with nothing in it and exit 0."""
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1000\t.\tA\tG\t50\tPASS\t.\n"
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode != 0
    assert "no sample columns" in result.stderr


def test_a_called_variant_without_an_ad_field_stops_the_run(tmp_path, fai):
    """VariantApplication indexes AlleleDepths[sample][alleleIndex] with no bounds check, and the
    array is empty when AD is absent - an IndexOutOfRangeException partway through a long run."""
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n"
        "1\t1000\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n"
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode != 0
    assert "no AD field" in result.stderr
    assert "bcftools +fill-tags" in result.stderr


def test_a_called_variant_with_a_dot_ad_stops_the_run(tmp_path, fai):
    """AD "." parses to a one-element array, so indexing it by allele 1 is out of range."""
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n"
        "1\t1000\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:.\n"
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode != 0
    assert "too short for allele 1" in result.stderr


def test_an_uncalled_sample_with_a_dot_ad_is_fine(tmp_path, fai):
    """What bcftools merge writes for a sample missing the variant. The builder skips it before
    reading any depth, so rejecting it would refuse every merged multi-sample VCF."""
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B\n"
        "1\t1000\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:10,8\t./.:.\n"
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode == 0, result.stderr
    assert "for 2 sample(s): SAMPLE_A, SAMPLE_B" in result.stderr


def test_multiallelic_padding_from_a_merge_is_fine(tmp_path, fai):
    """bcftools merge pads AD to Number=R with "." for alleles a sample does not carry. "." is a
    valid AD token to the builder and the array is still long enough to index."""
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B\n"
        "1\t1000\t.\tA\tG,T\t50\tPASS\t.\tGT:AD\t0/1:10,8,.\t0/2:12,.,6\n"
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode == 0, result.stderr


def test_the_reported_ad_problems_are_capped(tmp_path, fai):
    """A whole VCF from a caller that omits AD would otherwise print one line per variant."""
    lines = "".join(
        f"1\t{1000 + i}\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n" for i in range(40)
    )
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n" + lines
    )
    result = stage(tmp_path, fai, text)
    assert result.returncode != 0
    assert result.stderr.count("no AD field") <= 5
