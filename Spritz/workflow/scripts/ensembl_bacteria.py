"""Resolves and downloads an Ensembl Bacteria reference: genome FASTA, GFF3 and protein FASTA.

Bacteria are not on ftp.ensembl.org. They are served by Ensembl Genomes, on its own release
numbering - EG 63 corresponds to Ensembl 116 - and the layout differs from the vertebrate one in
three ways that between them rule out the URL templating downloads.smk uses for vertebrates:

1. There is an extra path level. A species lives under a numbered *collection*,
   `.../gff3/bacteria_5_collection/pseudomonas_aeruginosa_pao1_gca_000006765/`, and which collection
   is not derivable from the name: the 57 Pseudomonas aeruginosa strains are spread over about 30 of
   them. It has to be read from species_EnsemblBacteria.txt.

2. The DNA filename is not predictable. Most carry a trailing underscore on the assembly token -
   `...ASM676v1_.dna.toplevel.fa.gz` - and some do not, with nothing in the metadata to say which.
   So the filename is read from the small per-directory CHECKSUMS index rather than constructed.
   That also sidesteps assembly names needing sanitisation, which they otherwise do: 2,217 of the
   31,332 contain a space and 1,580 contain a '#'.

3. There is no `.dna.primary_assembly` dump, only `.dna.toplevel`.

The local filenames this writes are the ones the rest of the workflow already expects, including the
`.dna.primary_assembly.fa` name for what is really a toplevel dump. That is not a new inconsistency:
download_ensembl_references falls back from primary_assembly to toplevel for vertebrates too and
writes it under the same name.
"""

import argparse
import gzip
import re
import shutil
import sys
import urllib.request

FTP_ROOT = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/bacteria"
SPECIES_FILE = "species_EnsemblBacteria.txt"

# species_EnsemblBacteria.txt columns, 0-based, from its own '#name species division ...' header.
NAME = 0
SPECIES = 1
TAXONOMY_ID = 3
ASSEMBLY = 4
CORE_DB = 13

# core_db is e.g. "bacteria_5_collection_core_63_116_1"; the collection is everything before _core_.
COLLECTION = re.compile(r"^(.*)_core_\d+_\d+_\d+$")

# Ensembl derives a filename token from an assembly name by collapsing every run of characters
# outside this set to a single underscore, leaving dots and hyphens alone. Verified against the
# awkward cases in the metadata: "17870_2#15" -> "17870_2_15", "BRSU_AN4859/03" -> "BRSU_AN4859_03",
# "EC_O111:H8_CVM9634_1.0" -> "EC_O111_H8_CVM9634_1.0", "SBR5(T)" -> "SBR5_T_".
UNSAFE = re.compile(r"[^A-Za-z0-9.-]+")


def sanitise_assembly(assembly):
    """The filename-safe form of an assembly name, as Ensembl itself writes it."""
    return UNSAFE.sub("_", assembly)


def fetch(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read()


def release_root(release):
    return f"{FTP_ROOT}/release-{release}"


def species_file_url(release):
    return f"{release_root(release)}/{SPECIES_FILE}"


def parse_metadata(text):
    """species directory name -> (assembly, collection), for every row that has both."""
    rows = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) <= CORE_DB:
            continue
        collection = COLLECTION.match(fields[CORE_DB].strip())
        if not collection:
            continue
        rows[fields[SPECIES].strip()] = (fields[ASSEMBLY].strip(), collection.group(1))
    return rows


def parse_display_names(text):
    """species directory name -> the human-readable name in column 1."""
    names = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > SPECIES:
            names[fields[SPECIES].strip()] = fields[NAME].strip()
    return names


def parse_taxonomy_ids(text):
    """species directory name -> NCBI taxonomy id.

    Read by get_proteome.py: a bacterial proteome has to be looked up on UniProt by taxonomy id,
    because the name it is filed under is the strain-designation form rather than anything derivable
    from the Ensembl species name.
    """
    ids = {}
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > TAXONOMY_ID:
            ids[fields[SPECIES].strip()] = fields[TAXONOMY_ID].strip()
    return ids


def latest_release():
    """The newest Ensembl Genomes bacteria release, which is not the Ensembl release number.

    EG 63 is Ensembl 116. Only the EG numbering resolves on this FTP tree.
    """
    listing = fetch(f"{FTP_ROOT}/").decode("utf-8", errors="replace")
    found = {int(n) for n in re.findall(r"release-(\d+)", listing)}
    if not found:
        raise SystemExit(f"Error: no release-N directories found under {FTP_ROOT}/.")
    return str(max(found))


def resolve(metadata, species):
    """(assembly, collection) for a species directory name, or a SystemExit naming near misses."""
    found = metadata.get(species)
    if found:
        return found

    # The 2020-era bare names - pseudomonas_aeruginosa - no longer exist: every species directory now
    # carries a _gca_XXXXXXXXX accession suffix. So a miss is usually a name that used to be right,
    # and the useful reply is the strains that do exist.
    near = sorted(name for name in metadata if name.startswith(species))
    if not near:
        near = sorted(name for name in metadata if species in name)
    hint = (
        "Did you mean one of: " + ", ".join(near[:10]) + ("..." if len(near) > 10 else "")
        if near
        else "No species directory contains that string."
    )
    raise SystemExit(
        f"Error: Ensembl Bacteria has no species directory named '{species}'. Every species "
        f"directory carries a _gca_ accession suffix, so a bare species name will not match. {hint}"
    )


def checksum_names(url):
    """Filenames listed in a directory's CHECKSUMS index, whose third field is the name."""
    names = []
    for line in fetch(url).decode("utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 3:
            names.append(fields[-1])
    return names


def pick(names, suffix, where):
    """The one filename ending in `suffix`, erroring rather than guessing when that is not one."""
    matching = [name for name in names if name.endswith(suffix)]
    if len(matching) != 1:
        raise SystemExit(
            f"Error: expected exactly one file ending in '{suffix}' in {where}, found "
            f"{len(matching)}: {', '.join(sorted(matching)) or '(none)'}."
        )
    return matching[0]


def download_gunzipped(url, destination):
    with urllib.request.urlopen(url, timeout=600) as response:
        with gzip.GzipFile(fileobj=response) as unzipped:
            with open(destination, "wb") as out:
                shutil.copyfileobj(unzipped, out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, help="Ensembl Bacteria species directory name")
    parser.add_argument("--release", required=True, help="Ensembl Genomes release, e.g. 63")
    parser.add_argument("--assembly", help="expected assembly name, checked against the metadata")
    parser.add_argument("--metadata", help=f"a downloaded {SPECIES_FILE}; fetched if omitted")
    parser.add_argument("--genome-out", required=True)
    parser.add_argument("--gff3-out", required=True)
    parser.add_argument("--protein-out", required=True)
    args = parser.parse_args(argv)

    if args.metadata:
        with open(args.metadata, encoding="utf-8", errors="replace") as handle:
            text = handle.read()
    else:
        text = fetch(species_file_url(args.release)).decode("utf-8", errors="replace")

    species = args.species.lower()
    assembly, collection = resolve(parse_metadata(text), species)
    safe_assembly = sanitise_assembly(assembly)

    if args.assembly and args.assembly != safe_assembly:
        raise SystemExit(
            f"Error: the reference names assembly '{args.assembly}', but Ensembl Bacteria release "
            f"{args.release} has '{safe_assembly}' for {species}. Reference strings for bacteria are "
            f"release-specific; regenerate the row rather than editing it."
        )

    root = release_root(args.release)
    dna_dir = f"{root}/fasta/{collection}/{species}/dna"
    pep_dir = f"{root}/fasta/{collection}/{species}/pep"
    gff3_dir = f"{root}/gff3/{collection}/{species}"

    # Read the literal names rather than templating them. dna_rm./dna_sm. and .nonchromosomal are
    # excluded by the suffix; so are the .chromosome. GFF3 and the .abinitio. protein FASTA that
    # chromosome-level assemblies add alongside the ones wanted here.
    targets = [
        (f"{dna_dir}/{pick(checksum_names(dna_dir + '/CHECKSUMS'), '.dna.toplevel.fa.gz', dna_dir)}",
         args.genome_out),
        (f"{gff3_dir}/{pick(checksum_names(gff3_dir + '/CHECKSUMS'), f'.{args.release}.gff3.gz', gff3_dir)}",
         args.gff3_out),
        (f"{pep_dir}/{pick(checksum_names(pep_dir + '/CHECKSUMS'), '.pep.all.fa.gz', pep_dir)}",
         args.protein_out),
    ]

    for url, destination in targets:
        print(f"Downloading {url} -> {destination}", file=sys.stderr)
        download_gunzipped(url, destination)

    return 0


if __name__ == "__main__":
    sys.exit(main())
