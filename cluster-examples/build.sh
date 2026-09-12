#!/usr/bin/env bash
# Build the SIF on the cluster from the branch. Must run on x86_64: bioconda ships gatk4 for
# linux-64 only, so the variants environment will not solve on ARM.
set -euo pipefail
cd "$(dirname "$0")"

: "${SPRITZ_COMMIT:=vcf-entry-point}"
echo "Building spritz.sif from ${SPRITZ_COMMIT} on $(uname -m)"
[ "$(uname -m)" = "x86_64" ] || { echo "ERROR: must build on x86_64, this is $(uname -m)" >&2; exit 1; }

# --fakeroot avoids needing root. If the cluster disallows it, use a remote builder or ask an admin;
# the alternative (--sandbox without fakeroot) will not install conda packages correctly.
# The %test section in the definition runs at the end of the build, so a container that cannot
# find conda, snakemake or dotnet fails here instead of in an array task an hour later.
SPRITZ_COMMIT="$SPRITZ_COMMIT" apptainer build --fakeroot spritz.sif spritz.def

echo
echo "Built. It reports:"
apptainer exec spritz.sif cat /app/spritz/COMMIT
apptainer run --pwd /app/spritz/ spritz.sif \
  conda run --no-capture-output dotnet SpritzCMD.dll --help 2>&1 | head -5
