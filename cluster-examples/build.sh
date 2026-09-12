#!/usr/bin/env bash
# Build spritz.sif on the cluster from a branch, in two steps.
#
#   1. Compile SpritzCMD in a plain `apptainer exec` of the .NET SDK image. Unprivileged, writes to
#      this directory, and reproducible.
#   2. `apptainer build --fakeroot` a SIF that just copies that output in and installs the conda base
#      environment.
#
# The compile is deliberately NOT in the definition's %post. It was, and it failed under fakeroot
# with CS2001 on a file the SDK generates during build - which did not reproduce outside fakeroot.
# Moving it out removes the interaction rather than guessing at it, and makes step 2 fast to retry.
#
# x86_64 only: bioconda publishes gatk4 for linux-64 with no linux-aarch64, so the variants
# environment will not solve on ARM.
set -euo pipefail
cd "$(dirname "$0")"

: "${SPRITZ_COMMIT:=vcf-entry-point}"
: "${SDK_IMAGE:=docker://mcr.microsoft.com/dotnet/sdk:10.0}"

[ "$(uname -m)" = "x86_64" ] || { echo "ERROR: must build on x86_64, this is $(uname -m)" >&2; exit 1; }
command -v apptainer >/dev/null || { echo "ERROR: apptainer not found (module load apptainer?)" >&2; exit 1; }

# Somewhere with quota, per docs/wiki/Running-Spritz-on-a-cluster-with-Apptainer.md.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$PWD/.apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$PWD/.apptainer-tmp}"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "==> 1/3  fetching source at ${SPRITZ_COMMIT}"
rm -rf src build
git clone --quiet https://github.com/smith-chem-wisc/Spritz.git src
git -C src checkout --quiet "$SPRITZ_COMMIT"
git -C src rev-parse HEAD | tee src/COMMIT

echo "==> 2/3  compiling SpritzCMD in the .NET SDK container"
# HOME and the NuGet/CLI directories are redirected into this tree so nothing is written to a home
# directory that may be read-only or quota'd on a login node.
mkdir -p .dotnet-home .nuget
apptainer exec --cleanenv \
  --env HOME="$PWD/.dotnet-home" \
  --env DOTNET_CLI_HOME="$PWD/.dotnet-home" \
  --env NUGET_PACKAGES="$PWD/.nuget" \
  --env DOTNET_CLI_TELEMETRY_OPTOUT=1 \
  --env DOTNET_NOLOGO=1 \
  --bind "$PWD":"$PWD" --pwd "$PWD/src" \
  "$SDK_IMAGE" \
  dotnet build -c Release -p:UseSharedCompilation=false Spritz/SpritzCMD/SpritzCMD.csproj

OUT=src/Spritz/SpritzCMD/bin/Release/net10.0
test -f "$OUT/SpritzCMD.dll"        || { echo "ERROR: SpritzCMD.dll not produced" >&2; exit 1; }
test -f "$OUT/workflow/Snakefile"   || { echo "ERROR: workflow not copied into the build output" >&2; exit 1; }
mkdir -p build/app/spritz
cp -a "$OUT"/. build/app/spritz/
cp src/COMMIT build/app/spritz/COMMIT

echo "==> 3/3  building the SIF (the definition's %test runs at the end)"
apptainer build --fakeroot --force spritz.sif spritz.def

echo
echo "Built spritz.sif from $(cat build/app/spritz/COMMIT)"
