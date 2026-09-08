"""Tests for scripts/convert_ucsc2ensembl.py.

Renames UCSC-style chromosomes to Ensembl style and sorts the records karyotypically. Reads a
mapping file and config/config.yaml by relative path, reads a VCF on stdin and writes to stdout.
Invoked by rule download_dbsnp_vcf as:

    wget -O - <dbsnp> | zcat - | python scripts/convert_ucsc2ensembl.py > {output}
"""
from conftest import run_script

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"

MAPPING = "chr1\t1\nchr2\t2\nchr10\t10\nchrX\tX\nchrY\tY\nchrM\tMT\n"


def write_mapping(workflow_dir, contents=MAPPING, genome="GRCh38"):
    mapping = workflow_dir.parent / "resources" / "ChromosomeMappings" / f"{genome}_UCSC2ensembl.txt"
    mapping.write_text(contents)
    return mapping


def convert(workflow_dir, stdin):
    write_mapping(workflow_dir)
    result = run_script("convert_ucsc2ensembl.py", workflow_dir, stdin=stdin)
    assert result.returncode == 0, f"script failed: {result.stderr}"
    return result.stdout


def body(stdout):
    """Records only, dropping the header lines."""
    return [line for line in stdout.splitlines() if not line.startswith("#")]


def test_header_lines_pass_through_first(workflow_dir):
    out = convert(workflow_dir, HEADER + "chr1\t100\trs1\tA\tG\t.\t.\t.\n")
    assert out.startswith(HEADER)


def test_ucsc_names_are_renamed_to_ensembl(workflow_dir):
    out = convert(workflow_dir, HEADER + "chr1\t100\trs1\tA\tG\t.\t.\t.\n")
    assert body(out) == ["1\t100\trs1\tA\tG\t.\t.\t."]


def test_mitochondrion_is_renamed_from_chrM_to_MT(workflow_dir):
    out = convert(workflow_dir, HEADER + "chrM\t50\trsM\tA\tG\t.\t.\t.\n")
    assert body(out)[0].split("\t")[0] == "MT"


def test_records_are_emitted_in_karyotypic_order(workflow_dir):
    # Deliberately supplied out of order, and including 10 so that plain lexicographic sorting
    # would put it before 2.
    stdin = HEADER + "".join(
        f"{chrom}\t100\trs\tA\tG\t.\t.\t.\n" for chrom in ["chrX", "chr10", "chr2", "chrM", "chr1"]
    )
    out = convert(workflow_dir, stdin)
    assert [line.split("\t")[0] for line in body(out)] == ["1", "2", "10", "X", "MT"]


def test_unmapped_contigs_are_kept_after_the_named_chromosomes(workflow_dir):
    stdin = HEADER + "chrUn_something\t1\trs\tA\tG\t.\t.\t.\nchr1\t100\trs\tA\tG\t.\t.\t.\n"
    out = convert(workflow_dir, stdin)
    chroms = [line.split("\t")[0] for line in body(out)]
    assert chroms[0] == "1", "named chromosomes should come first"
    assert "chrUn_something" in chroms, "an unmapped contig should be retained, not dropped"


def test_script_leaves_the_rules_output_file_alone(workflow_dir):
    """Regression test for #253.

    Rule download_dbsnp_vcf redirects this script's stdout into
    ../resources/ensembl/<species>.ensembl.vcf. The script used to also open that same path with
    mode "w" and never write through the handle, truncating the file its caller was writing. It
    survived only because the open happened before any output, leaving the redirect's offset at 0.

    The script must not touch that path itself - everything goes to stdout.
    """
    ensembl_vcf = workflow_dir.parent / "resources" / "ensembl" / "Homo_sapiens.ensembl.vcf"
    ensembl_vcf.write_text("this content should survive a well-behaved script\n")

    out = convert(workflow_dir, HEADER + "chr1\t100\trs1\tA\tG\t.\t.\t.\n")

    assert body(out) == ["1\t100\trs1\tA\tG\t.\t.\t."], "records still go to stdout"
    assert ensembl_vcf.read_text() == "this content should survive a well-behaved script\n", (
        "the script must not open or truncate the path its rule redirects into"
    )
