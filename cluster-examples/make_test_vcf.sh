#!/usr/bin/env bash
# A small VCF against PAO1's own contig, for the `vcf` case. Deterministic and offline.
#
# The contig MUST be `Chromosome`: that is what Ensembl calls PAO1's single sequence, and a VCF
# using anything else is precisely the silent failure stage_user_vcf exists to catch (SnpEff reports
# ERROR_CHROMOSOME_NOT_FOUND and exits 0). Every sample carries AD, because the database builder
# indexes allele depths by allele number and a called variant without them is a hard error.
set -euo pipefail
out="${1:?usage: make_test_vcf.sh <output.vcf>}"
cat > "$out" <<'VCF'
##fileformat=VCFv4.2
##source=spritz-cluster-example
##contig=<ID=Chromosome,length=6264404>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE_A	SAMPLE_B
Chromosome	100000	.	A	G	60	PASS	.	GT:AD	0/1:22,18	./.:.
Chromosome	250000	.	C	T	60	PASS	.	GT:AD	1/1:0,31	0/1:14,9
Chromosome	500000	.	G	A	60	PASS	.	GT:AD	0/1:20,20	1/1:1,28
Chromosome	1000000	.	T	C	60	PASS	.	GT:AD	./.:.	0/1:11,13
Chromosome	2500000	.	A	T	60	PASS	.	GT:AD	0/1:19,17	0/1:15,16
VCF
echo "wrote $out"
