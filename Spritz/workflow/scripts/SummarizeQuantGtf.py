"""Combine StringTie transcript GTFs into one transcript-by-sample TPM matrix.

Joins on transcript ID, for the same reason as SummarizeQuantTab.py: the previous implementation
stacked values positionally and assumed every sample emitted its transcripts in the same order.
"""
import os
import re
import sys

import pandas as pd

# GTF attributes are `key "value"` pairs. Matching them rather than splitting on '; ' and ' '
# tolerates values containing spaces, missing trailing semicolons and extra attributes - none of
# which StringTie's -e output produces today, but a --merge GTF carrying gene_name would.
ATTRIBUTE = re.compile(r'(\w+)\s+"([^"]*)"')


def sample_name(path):
    """The sample label for a quant file, e.g. SRR13737862.sra.transcript.quant_ref.gtf -> SRR13737862."""
    return os.path.basename(path).split(".")[0]


def read_transcript_tpms(path):
    """Reads one StringTie GTF as a transcript-indexed TPM series.

    Only `transcript` features carry TPM; exon rows are skipped. Duplicated transcript IDs are
    summed, matching SummarizeQuantTab.py, so a repeated ID cannot silently overwrite an earlier
    value or break the join.
    """
    transcript_ids = []
    tpms = []
    with open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            # Raise rather than skip. A silent skip would turn a StringTie output change into an
            # empty or partial matrix that no downstream step notices; the version is not pinned,
            # so what produces these files can move under us.
            missing = [name for name in ("transcript_id", "TPM") if name not in attributes]
            if missing:
                raise ValueError(
                    f"{os.path.basename(path)} line {line_number}: transcript row has no "
                    f"{' or '.join(missing)}. Attributes present: "
                    f"{', '.join(sorted(attributes)) or 'none'}. "
                    "If StringTie changed its attribute names, this script needs updating."
                )
            transcript_ids.append(attributes["transcript_id"])
            tpms.append(float(attributes["TPM"]))

    if not transcript_ids:
        raise ValueError(
            f"{os.path.basename(path)}: no 'transcript' features found. Either the file is empty "
            "or StringTie no longer labels these rows 'transcript'."
        )

    frame = pd.DataFrame({"transcript_id": transcript_ids, "TPM": tpms})
    return frame.groupby("transcript_id")["TPM"].sum()


def build_matrix(input_files):
    """Transcripts down the rows, samples across the columns, aligned on transcript ID."""
    tpms_by_sample = {sample_name(path): read_transcript_tpms(path) for path in input_files}
    return pd.concat(tpms_by_sample, axis=1).sort_index().sort_index(axis=1)


def main():
    if len(sys.argv) < 3:
        print("Usage: SummarizeQuantGtf.py <output_file> <input_file_1> <input_file_2> ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    for path in input_files:
        print(f"reading {os.path.basename(path)}")

    print(f"Saving to {output_file} ...")
    build_matrix(input_files).to_csv(output_file)


if __name__ == "__main__":
    main()
