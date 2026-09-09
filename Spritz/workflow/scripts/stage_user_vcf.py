"""Stages a user-supplied VCF as the variant set, in place of calling variants from reads.

Two jobs beyond a copy.

Decompression: a VCF handed over by a collaborator is usually bgzipped, and the rest of the
workflow passes this file to SnpEff and to `cp` as plain text.

Contig checking: SnpEff reports a variant on a contig it does not know as ERROR_CHROMOSOME_NOT_FOUND
and carries on, so a VCF using UCSC names (chr1) against an Ensembl reference (1) yields a database
with every variant dropped and an exit code of 0. That is the failure this script exists to make
loud. It is a hard error only when *no* contig matches, because a VCF naming scaffolds the primary
assembly omits is normal and should not stop a run.
"""

import argparse
import gzip
import io
import sys


def read_fai_contigs(path):
    """Contig names from a samtools .fai, which is one tab-separated line per contig."""
    contigs = set()
    with open(path) as handle:
        for line in handle:
            name = line.split("\t")[0].strip()
            if name:
                contigs.add(name)
    return contigs


def open_maybe_gzip(path):
    """Reads .gz by content rather than by extension, since a .vcf is sometimes gzipped anyway."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def stage(vcf_path, fai_path, out):
    """Copies the VCF to `out`, returning (matched, unmatched) variant counts per contig name."""
    reference_contigs = read_fai_contigs(fai_path)
    matched = {}
    unmatched = {}

    with open_maybe_gzip(vcf_path) as handle:
        for line in handle:
            out.write(line)
            if line.startswith("#"):
                continue
            fields = line.split("\t", 1)
            contig = fields[0].strip()
            if not contig:
                continue
            tally = matched if contig in reference_contigs else unmatched
            tally[contig] = tally.get(contig, 0) + 1

    return matched, unmatched


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="the user-supplied VCF, optionally gzipped")
    parser.add_argument("fai", help="samtools .fai for the reference genome")
    args = parser.parse_args(argv)

    matched, unmatched = stage(args.vcf, args.fai, sys.stdout)
    matched_variants = sum(matched.values())
    unmatched_variants = sum(unmatched.values())

    if unmatched:
        # Capped: a VCF against the wrong assembly has thousands of distinct unknown contigs.
        examples = ", ".join(sorted(unmatched)[:10])
        print(
            f"{unmatched_variants} variant(s) on {len(unmatched)} contig(s) absent from the "
            f"reference will be dropped by SnpEff: {examples}",
            file=sys.stderr,
        )

    if not matched_variants:
        raise SystemExit(
            f"Error: no variant in {args.vcf} is on a contig named in the reference genome, so "
            f"annotating it would produce a database with no variants at all. The reference names "
            f"contigs like "
            f"{', '.join(sorted(read_fai_contigs(args.fai))[:5]) or '(none - the .fai is empty)'}, "
            f"while the VCF names them like {', '.join(sorted(unmatched)[:5]) or '(no variants)'}. "
            f"Ensembl-style names are expected; scripts/convert_ucsc2ensembl.py converts UCSC ones."
        )

    print(
        f"Staged {matched_variants} variant(s) on {len(matched)} reference contig(s).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
