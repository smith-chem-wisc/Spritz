"""Tests for scripts/simplify_gff3.py.

Strips exon and UTR features from a GFF3, so that SnpEff builds a database from transcript-level
records only. Takes the input path as argv[1] and writes to stdout. Invoked as:

    python scripts/simplify_gff3.py {input} > {output}
"""
from conftest import run_script

GFF3 = """\
##gff-version 3
##sequence-region 1 1 248956422
1\tensembl\tgene\t100\t500\t.\t+\t.\tID=gene:G1
1\tensembl\ttranscript\t100\t500\t.\t+\t.\tID=transcript:T1
1\tensembl\texon\t100\t200\t.\t+\t.\tID=exon:E1
1\tensembl\tCDS\t120\t200\t.\t+\t0\tID=CDS:C1
1\tensembl\tfive_prime_UTR\t100\t119\t.\t+\t.\tID=utr:U1
1\tensembl\tthree_prime_UTR\t450\t500\t.\t+\t.\tID=utr:U2
"""


def simplify(workflow_dir, contents=GFF3):
    gff = workflow_dir / "in.gff3"
    gff.write_text(contents)
    result = run_script("simplify_gff3.py", workflow_dir, args=[str(gff)])
    assert result.returncode == 0, f"script failed: {result.stderr}"
    return result.stdout


def feature_types(stdout):
    return [line.split("\t")[2] for line in stdout.splitlines() if len(line.split("\t")) > 2]


def test_exon_features_are_removed(workflow_dir):
    assert "exon" not in feature_types(simplify(workflow_dir))


def test_utr_features_are_removed(workflow_dir):
    types = feature_types(simplify(workflow_dir))
    assert not any(t.endswith("UTR") for t in types), types


def test_gene_transcript_and_cds_are_kept(workflow_dir):
    assert feature_types(simplify(workflow_dir)) == ["gene", "transcript", "CDS"]


def test_comment_lines_are_preserved(workflow_dir):
    """Regression test for #254.

    The guard used to be `emptyOrComment = len(linesplit) < 3 or line.startswith("#")`, with only
    lines where that was false written out - so every comment went, including the "##gff-version 3"
    pragma a GFF3 is supposed to open with. Removing exon and UTR features is the job; discarding
    the header was incidental.
    """
    out = simplify(workflow_dir)
    assert out.startswith("##gff-version 3\n"), "the pragma must stay on the first line"
    assert "##sequence-region 1 1 248956422" in out


def test_an_embedded_fasta_section_is_copied_through_intact(workflow_dir):
    """Sequence lines are not features and must not be filtered.

    SnpEff treats a "##FASTA" line as a switch to sequence-reading (FASTA_DELIMITER, gating
    readExonSequences). Emitting the marker while dropping the sequence under it - which the
    generic <3-column rule would do - is worse than dropping both: SnpEff would then read no exon
    sequences at all instead of falling back to the separate genome FASTA.
    """
    contents = GFF3 + "##FASTA\n>1\nACGTACGTAC\nGGGGCCCCAA\n"
    out = simplify(workflow_dir, contents)
    assert out.endswith("##FASTA\n>1\nACGTACGTAC\nGGGGCCCCAA\n"), out


def test_features_are_still_filtered_before_a_fasta_section(workflow_dir):
    contents = GFF3 + "##FASTA\n>1\nACGT\n"
    assert feature_types(simplify(workflow_dir, contents)) == ["gene", "transcript", "CDS"]


def test_a_feature_whose_name_merely_contains_utr_is_kept(workflow_dir):
    # The test is endswith("UTR"), so this is retained. Guards against a future change to a
    # substring match, which would silently discard records like this one.
    contents = GFF3 + "1\tensembl\tUTR_adjacent_region\t600\t700\t.\t+\t.\tID=x\n"
    assert "UTR_adjacent_region" in feature_types(simplify(workflow_dir, contents))
