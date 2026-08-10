"""Combine StringTie gene abundance tables (-A output) into one gene-by-sample TPM matrix.

Joins on gene ID. The previous implementation stacked the values positionally with numpy and
assumed every sample listed its genes in the same order, which is not true of StringTie output:
the same gene set comes back in a different order per sample. See the notes on each helper.
"""
import os
import sys

import pandas as pd


def sample_name(path):
    """The sample label for a quant file, e.g. SRR13737862.sra.gene.quant_ref.tab -> SRR13737862."""
    return os.path.basename(path).split(".")[0]


def read_gene_tpms(path):
    """Reads one StringTie -A table as a gene-indexed TPM series.

    StringTie can report the same gene ID on more than one line - it does so in the U2OS data, for
    gene:ENSG00000242759 - so the duplicates are summed rather than left to collide in the join.
    """
    table = pd.read_csv(path, sep="\t", usecols=["Gene ID", "TPM"])
    return table.groupby("Gene ID")["TPM"].sum()


def build_matrix(input_files):
    """Genes down the rows, samples across the columns, aligned on gene ID.

    pandas aligns on the index, so a sample that orders its genes differently still lands in the
    right rows. Both axes are sorted so the output does not depend on the order StringTie happened
    to emit, nor on the order snakemake happened to pass the files.
    """
    tpms_by_sample = {sample_name(path): read_gene_tpms(path) for path in input_files}
    return pd.concat(tpms_by_sample, axis=1).sort_index().sort_index(axis=1)


def main():
    if len(sys.argv) < 3:
        print("Usage: SummarizeQuantTab.py <output_file> <input_file_1> <input_file_2> ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]

    build_matrix(input_files).to_csv(output_file)

    print(f"Saved summarized data to {output_file}")


if __name__ == "__main__":
    main()
