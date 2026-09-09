import sys
gff = sys.argv[1]
with open(gff) as handle:
    for line in handle:
        # An embedded FASTA section is sequence, not tab-separated features. Copy it through
        # verbatim: the marker without its sequence would leave SnpEff reading no exon sequences.
        if line.startswith("##FASTA"):
            sys.stdout.write(line)
            sys.stdout.writelines(handle)
            break
        if line.startswith("#"):
            sys.stdout.write(line)
            continue
        linesplit = line.split('\t')
        if len(linesplit) < 3: continue
        if linesplit[2] == "exon" or linesplit[2].endswith("UTR"): continue
        sys.stdout.write(line)
