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
