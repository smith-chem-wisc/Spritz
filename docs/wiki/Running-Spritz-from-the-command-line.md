# Running Spritz from the command line

The recommended way to run Spritz without the GUI is to run the published container and let
`SpritzCMD` drive snakemake inside it. It needs a container runtime and nothing else — no conda, no
Python, no .NET, and no editing `config.yaml` by hand.

This is also the path with the strongest guarantee behind it: it is essentially what the
`dockerbuild` CI job runs on every push, so it is exercised continuously. The host-snakemake route at
the bottom of this page is not covered by any test.

## The short version

```bash
podman run --rm -it \
  -v "/path/to/analysis:/app/spritz/results/" \
  -v "/path/to/resources:/app/spritz/resources" \
  smithlab/spritz:0.3.14 \
  conda run --no-capture-output --live-stream \
  dotnet SpritzCMD.dll \
    -a=/app/spritz/results/ \
    -r="release-116,homo_sapiens,human,GRCh38" \
    -s=SRR629563 \
    -b -c -d
```

Substitute `docker` for `podman` if that is what you have; the flags are identical. On a cluster, see
[Running Spritz on a cluster with Apptainer](Running-Spritz-on-a-cluster-with-Apptainer).

## What the two mounts are for

Spritz writes into two directories and expects both to be supplied from the host, so results survive
the container exiting.

| Host | Inside the container | Holds |
|---|---|---|
| your analysis directory | `/app/spritz/results/` | FASTQs, alignments, the databases Spritz produces |
| your resources directory | `/app/spritz/resources` | Ensembl references, UniProt XML, the SnpEff database |

Point `resources` at somewhere you are willing to leave several tens of GB, and reuse the same one
across runs — the references are large and re-downloading them is the slowest part of a fresh run.

Note `-a=` takes the path **inside** the container (`/app/spritz/results/`), not your host path. The
mount is what connects them.

## The arguments

| Flag | Meaning |
|---|---|
| `-a=` | analysis directory, always `/app/spritz/results/` when running this way |
| `-r=` | reference, a line from `genomes.csv` — see below |
| `-s=` | paired-end SRA accession(s), comma-separated |
| `-t=` | single-end SRA accession(s) |
| `-f=` / `-i=` / `-j=` | local FASTQs instead of SRAs: single-end, first mate, second mate |
| `-b` | analyze variants |
| `-c` | analyze isoforms |
| `-d` | quantify |
| `-p=` | threads, defaults to the processor count |
| `--container-runtime` | `podman` (default), `docker`, or `apptainer` — only relevant when Spritz launches the container for you, not when you launch it yourself as above |

### Getting a reference string

The `-r=` value must be a line from `genomes.csv`, quoted, with four comma-separated fields:
`release,species,common name,assembly`. To see what is available:

```bash
podman run --rm -v "/path/to/analysis:/app/spritz/results/" smithlab/spritz:0.3.14 \
  conda run --no-capture-output dotnet SpritzCMD.dll -x -a=/app/spritz/results/
```

That writes `genomes.csv` into your analysis directory. Copy a line from it verbatim.

## Running a local build

To test a change to the workflow or to `SpritzCMD`, build the image yourself and run it by name:

```bash
podman build -t spritz:dev ./Spritz/
podman run --rm -it -v ... spritz:dev conda run ... dotnet SpritzCMD.dll ...
```

An image name that is not on a published registry is **not** pulled, so a local build is never
overwritten by the published one. If you publish to your own registry on a fork and *do* want it
pulled, that is what `RunnerEngine.AlwaysPull` is for.

## Where the results are

Everything lands under your analysis directory. The file to hand to MetaMorpheus is:

```
final/combined.spritz.snpeff.protein.withmods.xml.gz
```

with `final/combined.spritz.isoformvariants.protein.withmods.xml.gz` if you asked for isoforms as
well as variants.

## Advanced: running snakemake on the host

You can skip the container and run the workflow directly. This is the older documented route and it
is **not** covered by CI, so treat it as advanced.

```bash
conda env create --name spritzbase --file Spritz/workflow/envs/spritzbase.yaml
conda activate spritzbase
cd Spritz/workflow
snakemake -j 24 --use-conda --resources mem_mb=100000
```

Two things to know:

- **Do not pass `--conda-frontend mamba`.** Older instructions include it. Snakemake 9 accepts the
  flag, prints "Ignoring the alternative conda frontend setting", and uses conda — which now solves
  through libmamba anyway, so the flag only produces a warning.
- **You have to write `config/config.yaml` yourself.** The container path generates it from the
  command-line arguments; running snakemake directly does not.

### Native Windows

Running snakemake directly on native Windows fails during DAG construction with
`MissingInputException in rule all`, before any job starts. This is a path-separator bug, not a
problem with your installation — it is issue
[#243](https://github.com/smith-chem-wisc/Spritz/issues/243), fixed by
[#268](https://github.com/smith-chem-wisc/Spritz/pull/268).

Until that is released, on Windows use the container route above, or WSL. The container route is
unaffected, because the workflow runs on Linux inside the container regardless of your host.
