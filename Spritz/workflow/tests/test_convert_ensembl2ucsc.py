"""Tests for scripts/convert_ensembl2ucsc.py.

Renames Ensembl-style chromosomes to UCSC style and sorts the records. Takes input and output
paths as argv. Invoked by rule convert2ucsc.

Note that convert2ucsc is currently unreachable: nothing consumes its output
({dir}/isoforms/combined_ucsc.gtf), and it does not appear in the DAG for any analysis
combination. These tests therefore pin the behaviour of a script that does not presently run, so
that the traps in it are visible to whoever wires it back up.
"""
from conftest import run_script

MAPPING = "1\tchr1\n2\tchr2\n10\tchr10\nX\tchrX\nY\tchrY\nMT\tchrM\n"


def convert(workflow_dir, contents, mapping=MAPPING, mapping_genome="GRCh38"):
    mapping_path = (
        workflow_dir.parent / "resources" / "ChromosomeMappings" / f"{mapping_genome}_ensembl2UCSC.txt"
    )
    mapping_path.write_text(mapping)

    src = workflow_dir / "in.gtf"
    src.write_text(contents)
    dst = workflow_dir / "out.gtf"

    result = run_script("convert_ensembl2ucsc.py", workflow_dir, args=[str(src), str(dst)])
    return result, dst


def test_ensembl_names_are_renamed_to_ucsc(workflow_dir):
    result, dst = convert(workflow_dir, "1\tsource\texon\t1\t2\t.\t+\t.\tgene_id \"g\"\n")
    assert result.returncode == 0, result.stderr
    assert dst.read_text().split("\t")[0] == "chr1"


def test_mitochondrion_is_renamed_from_MT_to_chrM(workflow_dir):
    result, dst = convert(workflow_dir, "MT\tsource\texon\t1\t2\t.\t+\t.\tgene_id \"g\"\n")
    assert result.returncode == 0, result.stderr
    assert dst.read_text().split("\t")[0] == "chrM"


def test_records_are_emitted_in_karyotypic_order(workflow_dir):
    contents = "".join(
        f"{chrom}\tsource\texon\t1\t2\t.\t+\t.\tgene_id \"g\"\n" for chrom in ["X", "10", "2", "MT", "1"]
    )
    result, dst = convert(workflow_dir, contents)
    assert result.returncode == 0, result.stderr
    chroms = [line.split("\t")[0] for line in dst.read_text().splitlines()]
    assert chroms == ["chr1", "chr2", "chr10", "chrX", "chrM"]


def test_rows_with_an_empty_field_in_the_first_seven_are_dropped(workflow_dir):
    contents = (
        "1\tsource\texon\t1\t2\t.\t\t.\tgene_id \"bad\"\n"
        "2\tsource\texon\t1\t2\t.\t+\t.\tgene_id \"good\"\n"
    )
    result, dst = convert(workflow_dir, contents)
    assert result.returncode == 0, result.stderr
    assert "bad" not in dst.read_text()
    assert "good" in dst.read_text()


def test_the_configured_genome_is_ignored_in_favour_of_a_hardcoded_GRCh38(workflow_dir):
    """Pins a trap, and is the reason this file exists.

    The script reads "../resources/ChromosomeMappings/GRCh38_ensembl2UCSC.txt" as a literal, so
    the `genome` value in config/config.yaml has no effect. Its sibling
    convert_ucsc2ensembl.py interpolates {version} from config, so the two disagree.

    Consequence if convert2ucsc is ever wired back into the DAG: a non-GRCh38 run - yeast
    R64-1-1, mouse GRCm39 - would silently use *human* chromosome mappings. Names that are not in
    the human map fall through to `otherchr` and are written out unconverted, so the output would
    look plausible while carrying the wrong chromosome naming for the configured genome.

    This test asserts the current behaviour: a config naming R64-1-1 still reads the GRCh38 file,
    and fails if only an R64-1-1 mapping is present.
    """
    config = workflow_dir / "config" / "config.yaml"
    config.write_text(config.read_text().replace('genome: "GRCh38"', 'genome: "R64-1-1"'))

    # Only a correctly-named-for-the-config mapping file exists. A script honouring config would
    # find it; this one looks for GRCh38 and fails.
    result, _ = convert(workflow_dir, "I\tsource\texon\t1\t2\t.\t+\t.\tgene_id \"g\"\n",
                        mapping="I\tchrI\n", mapping_genome="R64-1-1")

    assert result.returncode != 0, "expected failure: the script ignores config and wants GRCh38"
    assert "GRCh38_ensembl2UCSC.txt" in result.stderr, result.stderr
