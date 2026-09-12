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


def called_alt_indices(genotype):
    """ALT allele indices a GT names, e.g. "0/1" -> {1}, "1|2" -> {1, 2}, "./." -> set()."""
    indices = set()
    for allele in genotype.replace("|", "/").split("/"):
        allele = allele.strip()
        if allele.isdigit() and allele != "0":
            indices.add(int(allele))
    return indices


def allele_depth_length(field):
    """How many AD entries the database builder will see.

    It splits on ',' discarding empty entries, so this has to discard them too - "10,,8" is two
    entries there, not three - and it indexes that array by allele index without a bounds check.
    """
    return len([part for part in field.split(",") if part.strip()])


def check_sample(format_keys, sample_field):
    """Why this sample would crash the database builder on this line, or None.

    VariantApplication reads AlleleDepths[sample][alleleIndex] directly, and the array is empty
    when AD is absent and length one when AD is ".". A sample whose genotype names an ALT allele
    is past the guard that skips uncalled samples, so a short array is an IndexOutOfRangeException
    partway through a long run rather than a dropped variant.
    """
    values = sample_field.split(":")
    fields = dict(zip(format_keys, values))
    called = called_alt_indices(fields.get("GT", "."))
    if not called:
        return None  # not called here, so its depths are never read
    if "AD" not in fields:
        return "no AD field"
    needed = max(called)
    if allele_depth_length(fields["AD"]) <= needed:
        return f"AD {fields['AD']!r} is too short for allele {needed}"
    return None


def stage(vcf_path, fai_path, out):
    """Copies the VCF to `out`, returning per-contig counts, sample names and AD problems."""
    reference_contigs = read_fai_contigs(fai_path)
    matched = {}
    unmatched = {}
    samples = []
    depth_problems = []

    with open_maybe_gzip(vcf_path) as handle:
        for line in handle:
            out.write(line)
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                # The one header line that is not a meta line names the samples, from column 10.
                samples = fields[9:]
                continue
            contig = fields[0].strip()
            if not contig:
                continue
            tally = matched if contig in reference_contigs else unmatched
            tally[contig] = tally.get(contig, 0) + 1

            if len(fields) > 9 and len(depth_problems) < 5:
                format_keys = fields[8].split(":")
                for name, sample_field in zip(samples, fields[9:]):
                    problem = check_sample(format_keys, sample_field)
                    if problem:
                        depth_problems.append(f"{contig}:{fields[1]} sample {name}: {problem}")
                        break

    return matched, unmatched, samples, depth_problems


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("vcf", help="the user-supplied VCF, optionally gzipped")
    parser.add_argument("fai", help="samtools .fai for the reference genome")
    args = parser.parse_args(argv)

    matched, unmatched, samples, depth_problems = stage(args.vcf, args.fai, sys.stdout)
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

    # A sites-only VCF is the other way to get a database with nothing in it: the database builder
    # keeps only variants that carry genotypes, so with no sample columns every one is discarded.
    if not samples:
        raise SystemExit(
            f"Error: {args.vcf} has no sample columns, so every variant in it would be discarded - "
            f"the database is built from genotypes, and a sites-only VCF has none. Supply a VCF with "
            f"at least one sample."
        )

    if depth_problems:
        raise SystemExit(
            "Error: variants are called in samples that have no usable AD (allele depth), which the "
            "database builder reads per allele:\n  "
            + "\n  ".join(depth_problems)
            + "\nAD has to come from the reads, so this needs re-genotyping rather than a tag fix: "
            "`bcftools mpileup -a AD` piped into `bcftools call`, or re-run the original caller. "
            "GATK emits AD by default; several other callers do not."
        )

    print(
        f"Staged {matched_variants} variant(s) on {len(matched)} reference contig(s) "
        f"for {len(samples)} sample(s): {', '.join(samples[:5])}"
        f"{'...' if len(samples) > 5 else ''}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
