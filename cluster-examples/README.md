# Bacterial verification on a cluster

The bacterial support in #280 has **never been executed** — the DAG resolves and the references
download, but SnpEff has never built a database from an Ensembl Bacteria reference and annotated
against it. That cannot be fixed on a Mac: bioconda publishes `gatk4` for `linux-64` only, with no
`linux-aarch64`, so the `variants` environment will not solve on ARM at all.

This runs it, from the branch rather than a published image, so the fixes are checked before merge.

## Order

```bash
./build.sh                 # x86_64 only; builds spritz.sif from the branch
sbatch run.slurm           # four cases as an array
./collect.sh               # tarball of logs, configs, prose and protein FASTAs
```

`build.sh` is two steps: it compiles `SpritzCMD` in a plain `apptainer exec` of the .NET SDK image,
then builds a SIF that only copies that output in and installs the conda environment. The compile
used to live in the definition's `%post`, which runs under `--fakeroot`, and failed there with a
`CS2001` on a file the SDK generates during build — a failure that did not reproduce outside
fakeroot. Moving the compile out removes the interaction instead of guessing at it, and means
iterating on the definition no longer recompiles.

The definition's `%test` section runs at the end of the build, so a container that cannot find
`conda`, `snakemake` or `dotnet` fails there rather than an hour into an array task. The first
version of this definition did not have it, and shipped with a `%runscript` that bypassed the base
image's entrypoint — so every invocation died with `exec: conda: not found`. The base image does not
put `/opt/conda/bin` on `PATH`; it relies on that entrypoint to activate the environment, which
Apptainer does not run when a definition supplies its own runscript.

`build.sh` needs `apptainer build --fakeroot`. The definition brings its own .NET SDK and clones the
branch, so the cluster needs neither `dotnet` nor Docker. Pin a commit with
`SPRITZ_COMMIT=<sha> ./build.sh` — the default is the branch tip, which moves.

## The cases

| | Organism | What it is for |
|---|---|---|
| `easy` | *P. aeruginosa* PAO1 | the happy path, plus bootstrapped known sites from real reads |
| `medium` | *Acetanaerobacterium elongatum* | assembly name with **spaces** — the case that breaks URL templating |
| `hard` | *Mycoplasma hyopneumoniae* | **codon table 4**, not 11 |
| `vcf` | PAO1 + a generated VCF | the `-v` entry point and its guards, offline and quick |

`medium` and `hard` request no analysis, so they build the reference database only. That is enough:
the codon table and the proteome lookup are what they test, and it avoids an SRA download.

**Why `hard` is the interesting one.** Mollicutes use NCBI table 4, where `TGA` is tryptophan rather
than a stop. Spritz emitted table 11 for every bacterium until this branch, which truncates every
protein at its first `TGA`. `check.sh` catches it without any proteomics, by comparing the protein
lengths SnpEff produces against the `pep.all.fa` Ensembl ships for the same genome — if the table is
wrong the median length collapses.

**Why `medium` is worth running.** 2,217 of Ensembl Bacteria's 31,332 assemblies contain a space and
1,580 contain a `#`. Ensembl rewrites those into filenames by collapsing unsafe runs to `_`, which is
why filenames are read from each directory's `CHECKSUMS` index rather than constructed. This is the
one case where getting that wrong produces a wrong URL rather than an obvious error.

## What `check.sh` asserts

Each check corresponds to something that is currently a structural claim only: the reference
downloaded from Ensembl Genomes, the codon table written into `snpEff.config`, **protein lengths
agreeing with Ensembl's own proteome**, the UniProt proteome resolved by taxonomy id, the config
recording `division` and `known_sites`, `prose.txt` describing what actually ran rather than a
`Homo_sapiens` quant, and the final database existing.

It runs automatically at the end of each array task, and can be re-run alone:

```bash
./check.sh hard work/hard
```

## Filesystems

`APPTAINER_TMPDIR` must be on a filesystem that supports **user extended attributes**. Rootless
`apptainer build --fakeroot` encodes file ownership in `user.rootlesscontainers` xattrs while
unpacking layers, and most shared HPC and network filesystems cannot store them — the build dies
with `unpriv.lsetxattr: invalid argument`, which does not mention the directory that caused it.

`build.sh` therefore defaults it to node-local scratch (`$SLURM_TMPDIR`, else `$TMPDIR`, else
`/tmp`) and probes it before starting, so a bad location fails in a second with an explanation rather
than several minutes in. Override it if your site puts scratch elsewhere:

```bash
APPTAINER_TMPDIR=/tmp/$USER/apptainer ./build.sh
```

The cache is only OCI blobs and needs no xattrs, so it stays in this directory where there is quota.

The base image is pulled to `micromamba-base.sif` once and the definition bootstraps from that rather
than from `docker://`. Converting an OCI image to SIF tolerates a filesystem without xattrs — it
warns and carries on — while the rootless unpack inside `build` does not, so pulling first removes
that failure mode rather than only relocating it. It also makes a rebuild skip the fetch.

## The warnings you will not see

`build.sh` silences two families, here only, because this is a verification build and both are
pre-existing on `master`:

- **`NU1902`** — `OpenMcdf 2.3.1` has moderate-severity advisories. No project references it
  directly; it arrives transitively through `mzLib 1.0.586`. Fixing it means bumping mzLib or pinning
  the transitive version, which needs its own testing and does not belong in a cluster script.
  `NuGetAudit=false` hides the report, not the risk — worth a separate issue.
- **`MSB3246`** — "PE image does not have metadata" while resolving references, from native
  libraries in the dependency set being handed to the reference resolver. Benign and long-standing.

The repository's own builds and CI are untouched, so the warning baseline the branch was measured
against still holds.

## Known ways this can fail that are not bugs

- **No UniProt proteome for the organism.** `get_proteome.py` looks bacteria up by NCBI taxonomy id
  and exits non-zero if UniProt has nothing for that taxon. Likely for obscure strains; `medium` is
  the one at risk. That is a genuine limitation, not a regression — say so rather than filing it.
- **`--fakeroot` disabled.** Some sites forbid it. There is no good workaround; a remote builder or
  an admin-built image is the answer.
- **SRA download rates.** `easy` fetches ~11M paired reads (DRR763891, verified present in ENA).

## Adding MetaMorpheus

Deliberately not wired in. The database is the thing under test here, and a search would add a large
dependency and a long runtime without testing anything this branch changed. Once a case passes, the
`final/*.protein.withmods.xml.gz` it produces is a normal MetaMorpheus input — and searching PAO1
proteomics against it, versus against the UniProt reference, is the first real measurement of whether
any of this finds peptides that reference-only search does not.
