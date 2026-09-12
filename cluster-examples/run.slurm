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

# NOT `cd "$(dirname "$0")"`. sbatch spools a *copy* of this script, so under SLURM $0 is something
# like /var/spool/slurmd/job1234/slurm_script and that cd lands in the spool directory, where
# spritz.sif, cases/ and check.sh do not exist. SLURM_SUBMIT_DIR is the directory sbatch was run
# from; fall back to $0 so the script is still runnable directly for a quick test.
HERE="${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
[ -f "$HERE/run.slurm" ] || HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
[ -f run.slurm ] || { echo "ERROR: cannot locate cluster-examples/ (tried $HERE); submit from that directory" >&2; exit 1; }
mkdir -p logs

CASES=(easy medium hard vcf)
CASE_NAME="${CASES[${SLURM_ARRAY_TASK_ID:-0}]}"
# shellcheck source=/dev/null
source "cases/${CASE_NAME}.env"

ROOT="${SPRITZ_WORK:-$PWD/work}/${CASE_NAME}"
RESULTS="$ROOT/results"; RESOURCES="$ROOT/resources"; SNAKE="$ROOT/snakemake"
mkdir -p "$RESULTS" "$RESOURCES" "$SNAKE" "$ROOT/cache"

# Point the apptainer cache at somewhere with quota, per the cluster wiki page.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$ROOT/../apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$ROOT/../apptainer-tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "=== ${CASE_NAME}: ${SPECIES} ${GENOME} (EG release ${RELEASE}) ==="
echo "    results=$RESULTS"

# `tr`, not ${SPECIES,,}: that lowercasing expansion is bash 4+, and it makes the script
# unrunnable on any bash 3.2 (macOS) where the rest of it could otherwise be tested.
SPECIES_LC="$(printf %s "$SPECIES" | tr '[:upper:]' '[:lower:]')"

ARGS=(-a=/app/spritz/results/
      -e=bacteria
      -p="${SLURM_CPUS_PER_TASK:-8}"
      -r="release-${RELEASE},${SPECIES_LC},${ORGANISM},${GENOME}")

# A case with reads calls variants, bootstrapping its known sites since Ensembl Bacteria has none.
# A case with a VCF annotates that instead. A case with neither builds the reference database only,
# which is enough to verify the download, the SnpEff build and the codon table.
if [ -n "${SRA:-}" ]; then ARGS+=(-s="$SRA" -b); fi
if [ -n "${VCF:-}" ]; then ./make_test_vcf.sh "$RESULTS/$VCF"; ARGS+=(-v="$VCF" -b); fi

# --pwd because Apptainer ignores the image's WORKDIR, so SpritzCMD.dll is not found without it.
# --cleanenv because a stray PYTHONPATH or CONDA_PREFIX from a loaded module otherwise shadows the
# image's own environment, which is the classic works-on-one-cluster-not-another failure.
#
# .snakemake is bound to real disk rather than left in the memory-backed --writable-tmpfs: snakemake
# builds a conda environment per rule there and they are not small, so a long run fills the tmpfs and
# dies partway through. Binding it also means the environments persist between the four cases instead
# of being rebuilt four times.
#
# `dotnet` directly, not via `conda run`. The definition's %environment runs the environment's
# activation hooks, so dotnet is on PATH already; and `conda run` is the machinery that exited 0
# while running nothing once a host conda leaked into the build. One less thing between here and the
# program. snakemake is still found on PATH by SpritzCMD, which is all it needs.
# snakemake keeps a source cache under $HOME/.cache/snakemake. --cleanenv leaves HOME pointing at
# the real home directory, which is writable, but it is set explicitly here so the cache cannot end
# up somewhere read-only - which is exactly how the container's own %test failed.
apptainer run --cleanenv --writable-tmpfs --pwd /app/spritz/ \
  --env XDG_CACHE_HOME="$ROOT/cache" \
  --bind "$RESULTS":/app/spritz/results/ \
  --bind "$RESOURCES":/app/spritz/resources \
  --bind "$SNAKE":/app/spritz/workflow/.snakemake \
  spritz.sif \
  dotnet SpritzCMD.dll "${ARGS[@]}"

echo "=== ${CASE_NAME}: run finished, checking ==="
./check.sh "$CASE_NAME" "$ROOT"
