#!/usr/bin/env bash
# Pull back what is worth looking at, not the genomes. Small enough to scp.
set -euo pipefail
cd "$(dirname "$0")"
# Includes the time: date alone meant a second collection the same day silently overwrote the
# first, which is exactly when you want both - before a fix and after it.
OUT="spritz-bacterial-results-$(date +%Y%m%d-%H%M%S).tar.gz"
tar czf "$OUT" \
  --exclude='*.fa' --exclude='*.fastq' --exclude='*.bam' --exclude='*.bai' --exclude='*.gz.tbi' \
  logs/ \
  $(find work -maxdepth 4 \
      \( -name conda-pkgs -o -name tmp -o -name cache \) -prune -o \
      \( -name 'prose.txt' -o -name 'config.yaml' -o -name '*.log' \
      -o -name 'snpEff.config' -o -name '*.protein.fasta' -o -name '*.vardesc.tsv' \
      -o -name '*.withmods.xml.gz' -o -name 'bootstrap.knownsites.vcf' \) -print 2>/dev/null) \
  2>/dev/null || true
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
