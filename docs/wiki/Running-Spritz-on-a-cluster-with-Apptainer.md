# Running Spritz on a cluster with Apptainer

Clusters almost never allow Docker, because it needs a daemon running as root. They usually provide
[Apptainer](https://apptainer.org) (formerly Singularity) instead, which runs unprivileged as you.
Spritz can drive it, using the same published image as the desktop runtimes.

## Check what your cluster has

```bash
apptainer --version    # or: singularity --version
module avail apptainer # many sites put it behind a module
```

If neither exists, ask your administrators — Apptainer is not something you can usefully install
yourself, since it needs setuid or unprivileged user namespaces configured by the site.

## Get the image

Prefer GHCR. Docker Hub rate-limits anonymous pulls **per IP address**, and every user on a login
node shares the same one, so Docker Hub pulls fail unpredictably on a busy cluster.

```bash
apptainer pull spritz.sif docker://ghcr.io/smith-chem-wisc/spritz:0.3.14
```

If your compute nodes have no outbound network, download the `.sif` attached to the
[release](https://github.com/smith-chem-wisc/Spritz/releases) on the login node and copy it across.

Build on a filesystem with room to spare, and point the cache somewhere with quota:

```bash
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer-cache
export APPTAINER_TMPDIR=/scratch/$USER/apptainer-tmp
```

## Run it

```bash
apptainer run --cleanenv --writable-tmpfs \
  --bind /scratch/$USER/spritz-analysis:/app/spritz/results/ \
  --bind /scratch/$USER/spritz-resources:/app/spritz/resources \
  spritz.sif \
  conda run --no-capture-output --live-stream dotnet SpritzCMD.dll \
    -a=/app/spritz/results/ \
    -r="release-116,homo_sapiens,human,GRCh38" \
    -s=SRR629563
```

Or let Spritz build the command for you:

```bash
SpritzCMD --container-runtime apptainer -a <analysis dir> -r "<reference>" -s <SRA>
```

### Why `--writable-tmpfs`

The image is read-only under Apptainer, and the workflow writes inside it — snakemake creates a conda
environment per rule at run time. Without this the run fails partway through rather than at startup.

**This matters for long runs.** That tmpfs is memory-backed and bounded, and those conda environments
are not small. If you see "no space left on device" from inside the container, bind a real directory
over the snakemake state instead:

```bash
  --bind /scratch/$USER/spritz-snakemake:/app/spritz/.snakemake
```

That also makes the environments persist between runs rather than being rebuilt each time.

### Why `--cleanenv`

Apptainer forwards your shell environment into the container by default, so a stray `PYTHONPATH`,
`CONDA_PREFIX` or `LD_LIBRARY_PATH` from a loaded module can shadow what is inside the image.
`--cleanenv` stops that. It is the single most common cause of a container that works on one cluster
and not another.

## Submitting it as a job

Apptainer runs in the foreground, so wrap it the way you would any other command. A SLURM example:

```bash
#!/bin/bash
#SBATCH --job-name=spritz
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

export APPTAINER_CACHEDIR=/scratch/$USER/apptainer-cache

apptainer run --cleanenv --writable-tmpfs \
  --bind /scratch/$USER/spritz-analysis:/app/spritz/results/ \
  --bind /scratch/$USER/spritz-resources:/app/spritz/resources \
  --bind /scratch/$USER/spritz-snakemake:/app/spritz/.snakemake \
  /scratch/$USER/spritz.sif \
  conda run --no-capture-output --live-stream dotnet SpritzCMD.dll \
    -a=/app/spritz/results/ -p=16 -r="release-116,homo_sapiens,human,GRCh38" -s=SRR629563
```

Ask for memory generously. Alignment and assembly dominate, and 24 GB is the desktop recommendation
for a human run.

## Differences from Docker and Podman

| | Docker / Podman | Apptainer |
|---|---|---|
| Daemon | yes (Docker) / no (Podman) | no |
| Mounts | `-v host:container` | `--bind host:container` |
| Container identity | named, `stop`-able | none; it is just a process |
| Image | pulled into a local store | a single `.sif` file you own |
| Runs as | root inside by default | **you**, always |

That last row is the one that surprises people: files Spritz writes are owned by *you*, not root.
This is usually an improvement — no more root-owned outputs — but it means the container cannot
write anywhere you cannot.

## Troubleshooting

**`FATAL: while extracting spritz.sif: root filesystem extraction failed`** — usually no space in
`APPTAINER_TMPDIR`. Point it at scratch.

**`no space left on device` partway through** — the `--writable-tmpfs` overlay filled. Bind a real
directory over `/app/spritz/.snakemake` as above.

**Permission denied writing results** — check the bound directory exists and is writable by you
before the run; Apptainer will not create it for you the way Docker does.

**Conda environments rebuilt on every run** — expected unless you bind `/app/spritz/.snakemake` to a
persistent directory. See [issue #36 in the tracker] for making this the default.
