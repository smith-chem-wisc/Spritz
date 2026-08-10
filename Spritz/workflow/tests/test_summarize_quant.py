"""Tests for SummarizeQuantTab.py and SummarizeQuantGtf.py, which build the final TPM matrices.

Fixtures mirror the U2OS run: Ensembl-style identifiers, StringTie's column header and filename
pattern, gene orders that differ between samples, and a duplicated gene ID.
"""
import csv

import pytest

from conftest import run_script

TAB_HEADER = "Gene ID\tGene Name\tReference\tStrand\tStart\tEnd\tCoverage\tFPKM\tTPM"

GENE_A = "gene:ENSG00000223972"
GENE_B = "gene:ENSG00000227232"
GENE_C = "gene:ENSG00000278267"
DUPLICATED_GENE = "gene:ENSG00000242759"

TRANSCRIPT_A = "transcript:ENST00000450305"
TRANSCRIPT_B = "transcript:ENST00000456328"
TRANSCRIPT_C = "transcript:ENST00000488147"


def write_tab(tmp_path, sample, genes_to_tpms):
    """Writes a StringTie -A gene abundance table, preserving the given gene order."""
    path = tmp_path / f"{sample}.fq.gene.quant_ref.tab"
    lines = [TAB_HEADER]
    for index, (gene, tpm) in enumerate(genes_to_tpms, start=1):
        lines.append(
            f"{gene}\t-\t1\t+\t{1000 * index}\t{1000 * index + 500}\t1.0\t{tpm / 2:.6f}\t{tpm:.6f}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def write_gtf(tmp_path, sample, transcripts_to_tpms, extra_lines=()):
    """Writes a StringTie transcript GTF, preserving the given transcript order."""
    path = tmp_path / f"{sample}.fq.transcript.quant_ref.gtf"
    lines = ["# stringtie -e -B -G reference.gff"]
    for index, (transcript, tpm) in enumerate(transcripts_to_tpms, start=1):
        start, end = 1000 * index, 1000 * index + 500
        lines.append(
            f'1\thavana\ttranscript\t{start}\t{end}\t.\t+\t.\t'
            f'gene_id "gene:ENSG0000000000{index}"; transcript_id "{transcript}"; '
            f'cov "1.0"; FPKM "{tpm / 2:.6f}"; TPM "{tpm:.6f}";'
        )
    lines.extend(extra_lines)
    path.write_text("\n".join(lines) + "\n")
    return path


def read_matrix(path):
    """Reads the written CSV as {row_id: {sample: float}}, plus the sample column order."""
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    samples = rows[0][1:]
    matrix = {row[0]: {s: float(v) for s, v in zip(samples, row[1:])} for row in rows[1:]}
    return samples, matrix


class TestSummarizeQuantTab:
    def test_keeps_every_gene_including_the_first(self, tmp_path):
        quant = write_tab(tmp_path, "SRR13737862", [(GENE_A, 20.1), (GENE_B, 6.83), (GENE_C, 162.4)])
        out = tmp_path / "gene.tpms.csv"

        result = run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(quant)])
        assert result.returncode == 0, result.stderr

        _, matrix = read_matrix(out)
        assert set(matrix) == {GENE_A, GENE_B, GENE_C}
        assert matrix[GENE_A]["SRR13737862"] == pytest.approx(20.1)

    def test_summarizes_multiple_samples(self, tmp_path):
        """The quant rules pass every sample at once, so this is the normal path."""
        first = write_tab(tmp_path, "sample_one", [(GENE_A, 1.0), (GENE_B, 2.0)])
        second = write_tab(tmp_path, "sample_two", [(GENE_A, 10.0), (GENE_B, 20.0)])
        out = tmp_path / "gene.tpms.csv"

        result = run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(first), str(second)])
        assert result.returncode == 0, result.stderr

        samples, matrix = read_matrix(out)
        assert samples == ["sample_one", "sample_two"]
        assert matrix[GENE_A] == {"sample_one": pytest.approx(1.0), "sample_two": pytest.approx(10.0)}

    def test_aligns_values_by_gene_id_not_row_position(self, tmp_path):
        """StringTie emits the same genes in a different order per sample, so the second sample here
        deliberately reorders them."""
        first = write_tab(tmp_path, "sample_one", [(GENE_A, 1.0), (GENE_B, 2.0), (GENE_C, 3.0)])
        second = write_tab(tmp_path, "sample_two", [(GENE_C, 30.0), (GENE_A, 10.0), (GENE_B, 20.0)])
        out = tmp_path / "gene.tpms.csv"

        run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(first), str(second)])

        _, matrix = read_matrix(out)
        assert matrix[GENE_A] == {"sample_one": pytest.approx(1.0), "sample_two": pytest.approx(10.0)}
        assert matrix[GENE_B] == {"sample_one": pytest.approx(2.0), "sample_two": pytest.approx(20.0)}
        assert matrix[GENE_C] == {"sample_one": pytest.approx(3.0), "sample_two": pytest.approx(30.0)}

    def test_sums_duplicated_gene_ids(self, tmp_path):
        """StringTie reports gene:ENSG00000242759 on two lines in the U2OS tables."""
        quant = write_tab(
            tmp_path, "SRR13737862",
            [(GENE_A, 1.0), (DUPLICATED_GENE, 0.882933), (DUPLICATED_GENE, 0.0)],
        )
        out = tmp_path / "gene.tpms.csv"

        result = run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(quant)])
        assert result.returncode == 0, result.stderr

        _, matrix = read_matrix(out)
        assert matrix[DUPLICATED_GENE]["SRR13737862"] == pytest.approx(0.882933)
        assert len(matrix) == 2, "a duplicated gene must collapse to one row, not appear twice"

    def test_sample_name_comes_from_the_filename(self, tmp_path):
        quant = write_tab(tmp_path, "1_130130_AH03CYADXX_P262_190E_index14", [(GENE_A, 1.0)])
        out = tmp_path / "gene.tpms.csv"

        run_script("SummarizeQuantTab.py", tmp_path, args=[str(out), str(quant)])

        samples, _ = read_matrix(out)
        assert samples == ["1_130130_AH03CYADXX_P262_190E_index14"], (
            "the sample label is the filename up to the first dot, and must not be truncated"
        )

    def test_usage_error_without_enough_arguments(self, tmp_path):
        result = run_script("SummarizeQuantTab.py", tmp_path, args=["only-the-output"])
        assert result.returncode == 1
        assert "Usage" in result.stdout


class TestSummarizeQuantGtf:
    def test_keeps_every_transcript_including_the_first(self, tmp_path):
        gtf = write_gtf(tmp_path, "SRR13737862",
                        [(TRANSCRIPT_A, 20.1), (TRANSCRIPT_B, 6.83), (TRANSCRIPT_C, 162.4)])
        out = tmp_path / "transcript.tpms.csv"

        result = run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])
        assert result.returncode == 0, result.stderr

        _, matrix = read_matrix(out)
        assert set(matrix) == {TRANSCRIPT_A, TRANSCRIPT_B, TRANSCRIPT_C}
        assert matrix[TRANSCRIPT_A]["SRR13737862"] == pytest.approx(20.1)

    def test_summarizes_multiple_samples(self, tmp_path):
        first = write_gtf(tmp_path, "sample_one", [(TRANSCRIPT_A, 1.0), (TRANSCRIPT_B, 2.0)])
        second = write_gtf(tmp_path, "sample_two", [(TRANSCRIPT_A, 10.0), (TRANSCRIPT_B, 20.0)])
        out = tmp_path / "transcript.tpms.csv"

        result = run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(first), str(second)])
        assert result.returncode == 0, result.stderr

        samples, matrix = read_matrix(out)
        assert samples == ["sample_one", "sample_two"]
        assert matrix[TRANSCRIPT_B] == {"sample_one": pytest.approx(2.0), "sample_two": pytest.approx(20.0)}

    def test_aligns_values_by_transcript_id_not_row_position(self, tmp_path):
        first = write_gtf(tmp_path, "sample_one",
                          [(TRANSCRIPT_A, 1.0), (TRANSCRIPT_B, 2.0), (TRANSCRIPT_C, 3.0)])
        second = write_gtf(tmp_path, "sample_two",
                           [(TRANSCRIPT_C, 30.0), (TRANSCRIPT_A, 10.0), (TRANSCRIPT_B, 20.0)])
        out = tmp_path / "transcript.tpms.csv"

        run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(first), str(second)])

        _, matrix = read_matrix(out)
        assert matrix[TRANSCRIPT_A] == {"sample_one": pytest.approx(1.0), "sample_two": pytest.approx(10.0)}
        assert matrix[TRANSCRIPT_C] == {"sample_one": pytest.approx(3.0), "sample_two": pytest.approx(30.0)}

    def test_ignores_comments_and_non_transcript_features(self, tmp_path):
        """Only `transcript` rows carry TPM."""
        exon = (
            '1\thavana\texon\t1000\t1500\t.\t+\t.\t'
            f'gene_id "gene:ENSG000000000001"; transcript_id "{TRANSCRIPT_A}"; cov "1.0";'
        )
        gtf = write_gtf(tmp_path, "SRR13737862", [(TRANSCRIPT_A, 20.1)], extra_lines=[exon])
        out = tmp_path / "transcript.tpms.csv"

        run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])

        _, matrix = read_matrix(out)
        assert set(matrix) == {TRANSCRIPT_A}
        assert matrix[TRANSCRIPT_A]["SRR13737862"] == pytest.approx(20.1)

    def test_parses_attribute_values_containing_spaces(self, tmp_path):
        """A --merge GTF carries gene_name, whose value can contain spaces."""
        line = (
            '1\thavana\ttranscript\t1000\t1500\t.\t+\t.\t'
            f'gene_id "gene:ENSG000000000001"; transcript_id "{TRANSCRIPT_A}"; '
            'ref_gene_name "some gene with spaces"; cov "1.0"; FPKM "10.0"; TPM "20.100000";'
        )
        gtf = tmp_path / "SRR13737862.fq.transcript.quant_ref.gtf"
        gtf.write_text(line + "\n")
        out = tmp_path / "transcript.tpms.csv"

        result = run_script("SummarizeQuantGtf.py", tmp_path, args=[str(out), str(gtf)])
        assert result.returncode == 0, result.stderr

        _, matrix = read_matrix(out)
        assert matrix[TRANSCRIPT_A]["SRR13737862"] == pytest.approx(20.1)

    def test_usage_error_without_enough_arguments(self, tmp_path):
        result = run_script("SummarizeQuantGtf.py", tmp_path, args=["only-the-output"])
        assert result.returncode == 1
        assert "Usage" in result.stdout
