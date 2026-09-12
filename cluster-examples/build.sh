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

# The cache is just OCI blobs, so it can live with the rest of this tree where there is quota.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$PWD/.apptainer-cache}"
mkdir -p "$APPTAINER_CACHEDIR"

# TMPDIR is different, and putting it here was wrong. `apptainer build --fakeroot` unpacks layers
# rootlessly, which encodes file ownership in `user.rootlesscontainers` extended attributes - and a
# shared HPC filesystem usually does not support xattrs, which fails the build with
#   unpriv.lsetxattr: invalid argument
# So this defaults to node-local scratch. Override APPTAINER_TMPDIR if your site puts scratch
# somewhere else, but it must support user xattrs.
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/apptainer-$USER}"
mkdir -p "$APPTAINER_TMPDIR"

# Checked rather than assumed, because the failure it produces names lsetxattr and not the directory.
if ! python3 - "$APPTAINER_TMPDIR" <<'XATTR'
import os, sys, tempfile
# Exit 0 = usable or undeterminable, 1 = definitely unsupported. Only a filesystem that answers
# "no" should block the build; not being able to ask is not a reason to refuse.
if not hasattr(os, "setxattr"):
    print("    (cannot probe xattrs on this platform, continuing)", file=sys.stderr)
    sys.exit(0)
try:
    with tempfile.NamedTemporaryFile(dir=sys.argv[1]) as f:
        os.setxattr(f.name, "user.spritz-probe", b"1")
except OSError as e:
    print(f"    {e.__class__.__name__}: {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"    (xattr probe inconclusive: {e}, continuing)", file=sys.stderr)
    sys.exit(0)
XATTR
then
    cat >&2 <<MSG

ERROR: $APPTAINER_TMPDIR does not support user extended attributes, which
       'apptainer build --fakeroot' needs in order to unpack image layers rootlessly.

       Point APPTAINER_TMPDIR at node-local scratch and re-run, for example:
           APPTAINER_TMPDIR=/tmp/\$USER/apptainer ./build.sh

       Shared project and network filesystems usually cannot do this; local disk can.
MSG
    exit 1
fi
echo "    APPTAINER_TMPDIR=$APPTAINER_TMPDIR (xattrs OK)"

echo "==> 1/3  fetching source at ${SPRITZ_COMMIT}"
rm -rf src build
git clone --quiet https://github.com/smith-chem-wisc/Spritz.git src
git -C src checkout --quiet "$SPRITZ_COMMIT"
git -C src rev-parse HEAD | tee src/COMMIT

echo "==> 2/3  compiling SpritzCMD in the .NET SDK container"
# HOME and the NuGet/CLI directories are redirected into this tree so nothing is written to a home
# directory that may be read-only or quota'd on a login node.
mkdir -p .dotnet-home .nuget
# No --env HOME: Apptainer refuses to override it ("Overriding HOME environment variable with
# APPTAINERENV_HOME is not permitted"). DOTNET_CLI_HOME and NUGET_PACKAGES are the two that decide
# where the SDK writes, and they are honoured.
apptainer exec --cleanenv \
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
# Bootstrap from a local SIF rather than docker://. Converting an OCI image to SIF tolerates a
# filesystem that cannot store xattrs - step 2 did exactly that, with warnings - while the rootless
# layer unpack that `build` performs does not. Pulling first therefore removes that failure mode
# rather than only relocating it, and makes a rebuild skip the fetch.
if [ ! -f micromamba-base.sif ]; then
    echo "    pulling the base image once"
    apptainer pull micromamba-base.sif docker://mambaorg/micromamba:2.9.0-ubuntu24.04
fi
apptainer build --fakeroot --force spritz.sif spritz.def

echo
echo "Built spritz.sif from $(cat build/app/spritz/COMMIT)"
