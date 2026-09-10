#!/usr/bin/env bash
# Pull back what is worth looking at, not the genomes. Small enough to scp.
set -euo pipefail
cd "$(dirname "$0")"
OUT="spritz-bacterial-results-$(date +%Y%m%d).tar.gz"
tar czf "$OUT" \
  --exclude='*.fa' --exclude='*.fastq' --exclude='*.bam' --exclude='*.bai' --exclude='*.gz.tbi' \
  logs/ \
  $(find work -maxdepth 4 \( -name 'prose.txt' -o -name 'config.yaml' -o -name '*.log' \
      -o -name 'snpEff.config' -o -name '*.protein.fasta' -o -name '*.vardesc.tsv' \
      -o -name '*.withmods.xml.gz' -o -name 'bootstrap.knownsites.vcf' \) 2>/dev/null) \
  2>/dev/null || true
echo "wrote $OUT ($(du -h "$OUT" | cut -f1))"
