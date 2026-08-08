"""Tests for SummarizeQuantTab.py and SummarizeQuantGtf.py.

These two scripts build the final TPM matrices - final/gene_reference_quant.tpms.csv and
final/transcript_reference_quant.tpms.csv - so they are the last step between StringTie's output
and what a user reads. They had no coverage, which is how the np.row_stack call below survived
until the Python upgrade removed it: NumPy 2.0 dropped row_stack, and the failure only surfaced
34 jobs into the container smoke run.

Two pre-existing defects are captured here as xfail rather than asserted as correct, so that
fixing them turns the tests green loudly instead of silently. See the tracking issue.
"""
import csv

from conftest import run_script

GENE_TAB = (
    "Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM\n"
    "YAL001C\tTFC3\tI\t-\t151006\t147594\t3.7\t12.5\t20.1\n"
    "YAL002W\tVPS8\tI\t+\t143707\t147531\t1.2\t4.25\t6.83\n"
    "YAL003W\tEFB1\tI\t+\t142174\t142253\t9.9\t101.0\t162.4\n"
)

TRANSCRIPT_GTF = (
    "# stringtie output\n"
    'I\tStringTie\ttranscript\t147594\t151006\t1000\t-\t.\tgene_id "YAL001C"; '
    'transcript_id "YAL001C_mRNA"; cov "3.7"; FPKM "12.500000"; TPM "20.100000";\n'
    'I\tStringTie\texon\t147594\t151006\t1000\t-\t.\tgene_id "YAL001C"; '
    'transcript_id "YAL001C_mRNA"; cov "3.7";\n'
    'I\tStringTie\ttranscript\t143707\t147531\t1000\t+\t.\tgene_id "YAL002W"; '
    'transcript_id "YAL002W_mRNA"; cov "1.2"; FPKM "4.250000"; TPM "6.830000";\n'
    'I\tStringTie\ttranscript\t142174\t142253\t1000\t+\t.\tgene_id "YAL003W"; '
    'transcript_id "YAL003W_mRNA"; cov "9.9"; FPKM "101.000000"; TPM "162.400000";\n'
)


def read_matrix(path):
    """Reads the written CSV into (sample_columns, {row_id: {sample: value}})."""
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    samples = rows[0][1:]
    return samples, {row[0]: dict(zip(samples, row[1:])) for row in rows[1:]}


class TestSummarizeQuantTab:
    def test_runs_and_writes_a_matrix(self, tmp_path):
        """Regression guard for the NumPy 2.0 removal of np.row_stack.

        Before the fix this died with AttributeError: module 'numpy' has no attribute
        'row_stack', which is what took down make_gene_quant_dataframe_ref.
        """
        quant = tmp_path / "SRR13737862.sra.gene.quant_ref.tab"
        quant.write_text(GENE_TAB)
        out = tmp_path / "gene_reference_quant.tpms.csv"

        result = run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(quant)])

        assert result.returncode == 0, result.stderr
        assert out.exists()
        samples, matrix = read_matrix(out)
        assert samples == ["SRR13737862"]
        assert matrix["YAL003W"]["SRR13737862"] == "162.4"

    def test_usage_error_without_enough_arguments(self, tmp_path):
        result = run_script("SummarizeQuantTab.py", tmp_path, args=["only-the-output"])
        assert result.returncode == 1
        assert "Usage" in result.stdout

    def test_keeps_every_gene(self, tmp_path):
        """KNOWN BROKEN: the first gene of every file is dropped.

        `dataframe[1:, 0] = [sample names]` overwrites column 0 - which holds real data,
        because no placeholder cell was prepended to the ID row - and the subsequent
        `[1:, 1:]` slice then discards it. So the first gene's TPM is destroyed and the
        gene disappears from the output entirely.
        """
        quant = tmp_path / "SRR13737862.sra.gene.quant_ref.tab"
        quant.write_text(GENE_TAB)
        out = tmp_path / "gene_reference_quant.tpms.csv"

        run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(quant)])

        _, matrix = read_matrix(out)
        if "YAL001C" not in matrix:
            import pytest

            pytest.xfail("first gene dropped; see the SummarizeQuant data-loss issue")
        assert matrix["YAL001C"]["SRR13737862"] == "20.1"


class TestSummarizeQuantGtf:
    def test_runs_and_writes_a_matrix(self, tmp_path):
        """Regression guard for the NumPy 2.0 removal of np.row_stack."""
        gtf = tmp_path / "SRR13737862.sra.transcript.quant_ref.gtf"
        gtf.write_text(TRANSCRIPT_GTF)
        out = tmp_path / "transcript_reference_quant.tpms.csv"

        result = run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])

        assert result.returncode == 0, result.stderr
        assert out.exists()
        samples, matrix = read_matrix(out)
        assert samples == ["SRR13737862"]
        assert matrix["YAL003W_mRNA"]["SRR13737862"] == "162.400000"

    def test_ignores_non_transcript_rows(self, tmp_path):
        """Only `transcript` features carry TPM; exon rows must not become columns."""
        gtf = tmp_path / "SRR13737862.sra.transcript.quant_ref.gtf"
        gtf.write_text(TRANSCRIPT_GTF)
        out = tmp_path / "transcript_reference_quant.tpms.csv"

        run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])

        _, matrix = read_matrix(out)
        assert all(not key.endswith("exon") for key in matrix)
        assert len(matrix) <= 3

    def test_keeps_every_transcript(self, tmp_path):
        """KNOWN BROKEN: same first-record loss as the .tab script."""
        gtf = tmp_path / "SRR13737862.sra.transcript.quant_ref.gtf"
        gtf.write_text(TRANSCRIPT_GTF)
        out = tmp_path / "transcript_reference_quant.tpms.csv"

        run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])

        _, matrix = read_matrix(out)
        if "YAL001C_mRNA" not in matrix:
            import pytest

            pytest.xfail("first transcript dropped; see the SummarizeQuant data-loss issue")
        assert matrix["YAL001C_mRNA"]["SRR13737862"] == "20.100000"
