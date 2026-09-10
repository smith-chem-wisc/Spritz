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
  smithlab/spritz:0.3.15 \
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
| `-v=` | comma-separated VCFs you called elsewhere, annotated instead of calling variants from reads — see below |
| `-e=` | Ensembl division: `vertebrates` (default) or `bacteria` — see below |
| `-b` | analyze variants |
| `-c` | analyze isoforms |
| `-d` | quantify |
| `-p=` | threads, defaults to the processor count |
| `--container-runtime` | `podman` (default), `docker`, or `apptainer` — only relevant when Spritz launches the container for you, not when you launch it yourself as above |

### Annotating a VCF you already have

If variants were called outside Spritz — from WGS or exome reads, or by any caller — `-v=` skips
alignment and GATK entirely and annotates that VCF directly:

```bash
podman run --rm -it \
  -v "/path/to/analysis:/app/spritz/results/" \
  -v "/path/to/resources:/app/spritz/resources" \
  smithlab/spritz:0.3.15 \
  conda run --no-capture-output --live-stream \
  dotnet SpritzCMD.dll \
    -a=/app/spritz/results/ \
    -r="release-116,homo_sapiens,human,GRCh38" \
    -v=my_variants.vcf \
    -b
```

Three things to know.

**`-v=` takes filenames, not paths**, resolved inside your analysis directory — the same convention
as `-i=`/`-j=`/`-f=`. Only the analysis and resources directories are mounted into the container, so
a host path from anywhere else would not resolve. Copy or move the VCFs into the analysis directory
first. Gzipped VCFs are accepted.

**One VCF per sample is the expected case.** Pass them comma-separated —
`-v=sample_a.vcf,sample_b.vcf,sample_c.vcf` — and Spritz merges them into a single multi-sample VCF
with `bcftools merge` before annotating. That is worth doing rather than concatenating, because the
database builder is sample-aware: it walks the genotypes per individual and emits variant protein
sequences for each, so a merged multi-sample VCF gives you per-individual variant proteins. Note this
is *more* faithful than the read-based path, which assigns every input the same read group and pools
everything into one sample.

**Your VCFs must carry genotypes and allele depths.** Two requirements, both checked before
anything runs, because both otherwise fail late and unhelpfully:

- **At least one sample column.** The database is built from genotypes, so a sites-only VCF yields a
  database with nothing in it — and would previously have done so while exiting 0.
- **A per-sample `AD` (allele depth) for every variant a sample actually carries.** The builder
  indexes allele depths by allele number, so a called variant with no `AD` is an error partway
  through the run. GATK emits `AD` by default; several other callers do not. If yours does not, add
  it with `bcftools +fill-tags -- -t AD`. A sample that simply does not carry a variant needs
  nothing — `./.:.` is what a merge writes there and is fine.

Sample names that collide between files are renamed rather than rejected, since callers often emit a
placeholder name. The renaming is positional (`SAMPLE`, `2:SAMPLE`, `3:SAMPLE`), so if you want the
merged columns to identify their source, give each VCF a distinct sample name before passing it in.

With a single VCF none of this applies: the file goes straight to annotation untouched.

**It requires `-b` and excludes `-c` and `-d`.** Isoform reconstruction assembles transcripts and
quantification counts reads, so neither has an input without them. It also cannot be combined with
`-s=`/`-t=`/`-i=`/`-j=`/`-f=`: a supplied VCF replaces variant calling, so passing both is ambiguous
and is rejected rather than silently resolved.

**Contig names must match the Ensembl reference.** Ensembl calls the first human chromosome `1`; a
VCF from a UCSC-based pipeline calls it `chr1`. SnpEff reports a variant on a contig it does not know
as `ERROR_CHROMOSOME_NOT_FOUND` and exits 0, so a mismatch would otherwise hand you a database with
every variant silently dropped. Spritz checks this before annotating and stops the run when nothing
matches. A VCF that merely names some scaffolds the primary assembly omits is fine — those are
reported and skipped.

The reference database itself is still built from Ensembl, so the genome, GFF3 and protein FASTA are
downloaded as usual; what `-v=` saves is the read download, trimming, alignment and variant calling.

### Bacterial references

Bacteria are not on `ftp.ensembl.org`. They come from Ensembl Genomes, which numbers its releases
separately — **EG 63 is Ensembl 116** — so a bacterial run needs `-e=bacteria` and an EG release
number:

```bash
dotnet SpritzCMD.dll \
  -a=/app/spritz/results/ \
  -e=bacteria \
  -r="release-63,pseudomonas_aeruginosa_pao1_gca_000006765,pseudomonas aeruginosa pao1,ASM676v1" \
  -v=my_variants.vcf \
  -b
```

**A bacterial reference requires `-v=`.** Ensembl Bacteria publishes no known variant sites — its
`variation/` directory holds only a VEP cache, with no `vcf/` — and GATK base recalibration, the only
thing that reads them, has no input without them. So Spritz cannot call bacterial variants from reads,
and rejects the combination up front rather than failing partway through a long run. Call the variants
with your own pipeline and hand Spritz the VCF.

#### Finding the reference string

Species directory names are strain-specific and carry a GCA accession —
`pseudomonas_aeruginosa_pao1_gca_000006765`, not `pseudomonas_aeruginosa` — and there are 31,332 of
them, so `-x` does not list them. Search for yours:

```bash
podman run --rm -v "/path/to/analysis:/app/spritz/results/" smithlab/spritz:0.3.14 \
  conda run --no-capture-output python workflow/scripts/update_genomes.py \
    --division bacteria --match pseudomonas_aeruginosa \
    --output /app/spritz/results/genomes.csv
```

That appends the matching rows to `genomes.csv` in your analysis directory, alongside any vertebrate
rows already there. Copy one out verbatim. Passing a bare species name to `-r=` instead will fail with
an error listing the strains that do exist.

#### Two things that differ from a vertebrate run

**Codon table.** Bacterial genomes are translated with NCBI table 11
(`Bacterial_and_Plant_Plastid`) rather than the standard table. Against the standard table it differs
in exactly four codons — `ATT`, `ATC`, `ATA` and `GTG` — and in each only by whether the codon counts
as a valid start. Every codon-to-amino-acid mapping is the same, so this affects `start_lost` and
initiation calls and nothing else; missense, synonymous and stop_gained are identical either way.
Note that SnpEff ships its own bacterial genome entries and *none* of them declares a codon table, so
they all translate with the standard one; the database Spritz builds here does declare it.

**Contig names.** A bacterial assembly typically has one sequence, and Ensembl names it `Chromosome`,
not `1`. Your VCF's `CHROM` column has to match — see the contig-name note above.

### Getting a reference string

The `-r=` value must be a line from `genomes.csv`, quoted, with four comma-separated fields:
`release,species,common name,assembly`. To see what is available:

```bash
podman run --rm -v "/path/to/analysis:/app/spritz/results/" smithlab/spritz:0.3.15 \
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
