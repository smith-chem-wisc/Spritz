"""Tests for scripts/clean_vcf.py, which filters a VCF stream on stdin.

Pure stdin to stdout with no configuration or resource files, so it is the most directly
testable script in the workflow.
"""
from conftest import run_script

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def clean(stdin, cwd):
    result = run_script("clean_vcf.py", cwd, stdin=stdin)
    assert result.returncode == 0, f"clean_vcf.py failed: {result.stderr}"
    return result.stdout


def test_header_lines_pass_through_unchanged(workflow_dir):
    assert clean(HEADER, workflow_dir) == HEADER


def test_ordinary_variant_is_kept(workflow_dir):
    row = "1\t100\trs1\tA\tG\t50\tPASS\tAC=1\n"
    assert clean(HEADER + row, workflow_dir) == HEADER + row


def test_row_with_an_empty_field_in_the_first_seven_is_dropped(workflow_dir):
    # An empty ID column. The script rejects an empty value anywhere in columns 0-6.
    bad = "1\t100\t\tA\tG\t50\tPASS\tAC=1\n"
    good = "1\t200\trs2\tC\tT\t50\tPASS\tAC=1\n"
    assert clean(HEADER + bad + good, workflow_dir) == HEADER + good


def test_alt_allele_without_a_recognised_base_is_dropped(workflow_dir):
    # ALT of "<DEL>" shares no character with {A,C,T,G}, so it is unparsable here.
    symbolic = "1\t100\trs1\tA\t<DEL>\t50\tPASS\tSVTYPE=DEL\n"
    kept = "1\t200\trs2\tA\tG\t50\tPASS\tAC=1\n"
    assert clean(HEADER + symbolic + kept, workflow_dir) == HEADER + kept


def test_alt_with_more_than_one_duplicate_allele_pattern_is_dropped(workflow_dir):
    # The script drops rows whose ALT contains more than one "<letter>." pattern.
    duplicated = "1\t100\trs1\tA\tG.,T.\t50\tPASS\tAC=1\n"
    kept = "1\t200\trs2\tA\tG\t50\tPASS\tAC=1\n"
    assert clean(HEADER + duplicated + kept, workflow_dir) == HEADER + kept


def test_a_single_duplicate_allele_pattern_is_still_kept(workflow_dir):
    # Boundary of the rule above: one occurrence is allowed, two are not.
    single = "1\t100\trs1\tA\tG.\t50\tPASS\tAC=1\n"
    assert clean(HEADER + single, workflow_dir) == HEADER + single


def test_empty_input_produces_empty_output(workflow_dir):
    assert clean("", workflow_dir) == ""
