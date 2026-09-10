#!/usr/bin/env bash
#SBATCH --job-name=spritz-bact
#SBATCH --array=0-3
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%a.out
#SBATCH --error=logs/%x-%a.err
#
# Four bacterial cases as an array. Adjust the SBATCH lines for your partition/account.
#
#   sbatch run.slurm
#   sbatch --array=3 run.slurm      # just the hard (codon table) case
#
# GATK asks for -Xmx24000M in the calling rules, so 32G is the floor for the case that calls
# variants. The three reference-only cases are far lighter but share the allocation.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

CASES=(easy medium hard vcf)
CASE_NAME="${CASES[${SLURM_ARRAY_TASK_ID:-0}]}"
# shellcheck source=/dev/null
source "cases/${CASE_NAME}.env"

ROOT="${SPRITZ_WORK:-$PWD/work}/${CASE_NAME}"
RESULTS="$ROOT/results"; RESOURCES="$ROOT/resources"
mkdir -p "$RESULTS" "$RESOURCES" 

echo "=== ${CASE_NAME}: ${SPECIES} ${GENOME} (EG release ${RELEASE}) ==="
echo "    results=$RESULTS"

ARGS=(-a=/app/spritz/results/
      -e=bacteria
      -p="${SLURM_CPUS_PER_TASK:-8}"
      -r="release-${RELEASE},${SPECIES,,},${ORGANISM},${GENOME}")

# A case with reads calls variants, bootstrapping its known sites since Ensembl Bacteria has none.
# A case with a VCF annotates that instead. A case with neither builds the reference database only,
# which is enough to verify the download, the SnpEff build and the codon table.
if [ -n "${SRA:-}" ]; then ARGS+=(-s="$SRA" -b); fi
if [ -n "${VCF:-}" ]; then ./make_test_vcf.sh "$RESULTS/$VCF"; ARGS+=(-v="$VCF" -b); fi

apptainer run --cleanenv --writable-tmpfs \
  --bind "$RESULTS":/app/spritz/results/ \
  --bind "$RESOURCES":/app/spritz/resources \
  spritz.sif dotnet SpritzCMD.dll "${ARGS[@]}"

echo "=== ${CASE_NAME}: run finished, checking ==="
./check.sh "$CASE_NAME" "$ROOT"
