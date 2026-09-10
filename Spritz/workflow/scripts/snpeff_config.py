"""Builds the genome entry appended to SnpEff's snpEff.config before a database is built.

SnpEff translates every chromosome with the standard codon table unless the config names another
one for it. CodonTables.getTable(genome, chromosome) looks up genome+chromosome, then genome, then
falls back to Standard; there is no name-based detection of a mitochondrial contig. So whatever
this emits is the whole story for mitochondrial protein-coding genes, and getting it wrong changes
SnpEff's effect calls on them - synonymous against missense against stop_gained - which is what
feeds the variant protein XML.

Ensembl's main site is vertebrates, whose mitochondria use NCBI translation table 2, so that is the
default. The exceptions are the few non-vertebrates Ensembl carries. Each of those also names its
mitochondrial contig something other than MT, so a hardcoded MT line does not merely give them the
wrong table, it does not apply to them at all and leaves their mitochondrial genes on Standard,
where TGA is a stop rather than tryptophan.

Contig names verified against the Ensembl REST assembly info for each species.
"""

# Names defined by codon.* entries in the SnpEff config, so a table named here resolves rather than
# silently falling back. Not the full list SnpEff ships - only what this module can emit.
KNOWN_CODON_TABLES = frozenset(
    {
        "Ascidian_Mitochondrial",
        "Bacterial_and_Plant_Plastid",
        "Mycoplasma",
        "Spiroplasma",
        "Invertebrate_Mitochondrial",
        "Vertebrate_Mitochondrial",
        "Yeast_Mitochondrial",
    }
)

# NCBI translation table 11, declared for the whole genome rather than one contig: a bacterium has a
# single chromosome and no mitochondrion, and CodonTables.getTable falls back genome+chromosome ->
# genome -> Standard, so a genome-level entry covers every contig including plasmids.
#
# Against codon.Standard in the shipped snpEff.config it differs in exactly four codons, and in each
# only by the "+" that marks a valid start: ATT/I, ATC/I, ATA/I and GTG/V. The codon-to-amino-acid
# map is identical. So this changes start_lost and initiation calls and nothing else - missense,
# synonymous and stop_gained are unaffected either way.
BACTERIAL_CODON_TABLE = "Bacterial_and_Plant_Plastid"

# Not every bacterium uses table 11. The Mollicutes - Mycoplasma and its relatives, and Spiroplasma -
# use NCBI table 4, and the difference is not a nuance: against table 11 they differ at TGA, which is
# a stop under 11 and tryptophan under 4. Translating one of these genomes with table 11 truncates
# every protein at its first TGA.
#
# Matched on the genus prefix of the Ensembl species directory name, which is a heuristic rather
# than a taxonomy lookup: NCBI assigns table 4 across Mycoplasmatales and Entomoplasmatales, and the
# genera below are the ones Ensembl Bacteria actually carries. Getting it wrong is detectable rather
# than silent - the protein lengths SnpEff emits stop matching the pep.all.fa Ensembl ships for the
# same genome - which is what the bacterial verification case checks.
MYCOPLASMA_CODON_TABLE = "Mycoplasma"
SPIROPLASMA_CODON_TABLE = "Spiroplasma"
TABLE_4_GENERA = {
    "mycoplasma": MYCOPLASMA_CODON_TABLE,
    "mycoplasmoides": MYCOPLASMA_CODON_TABLE,
    "mycoplasmopsis": MYCOPLASMA_CODON_TABLE,
    "mesoplasma": MYCOPLASMA_CODON_TABLE,
    "entomoplasma": MYCOPLASMA_CODON_TABLE,
    "ureaplasma": MYCOPLASMA_CODON_TABLE,
    "malacoplasma": MYCOPLASMA_CODON_TABLE,
    "metamycoplasma": MYCOPLASMA_CODON_TABLE,
    "spiroplasma": SPIROPLASMA_CODON_TABLE,
}


def bacterial_codon_table(species):
    """The codon table for a bacterial genome, which is not always table 11."""
    genus = species.lower().split("_")[0]
    return TABLE_4_GENERA.get(genus, BACTERIAL_CODON_TABLE)

# Ensembl serves bacteria from Ensembl Genomes, on its own release numbering, not from ftp.ensembl.org.
ENSEMBL_REFERENCE = "https://ftp.ensembl.org/pub/release-{release}/"
ENSEMBL_BACTERIA_REFERENCE = "https://ftp.ebi.ac.uk/ensemblgenomes/pub/bacteria/release-{release}/"

# M is inert for Ensembl, which spells the contig MT. Kept because dropping it changes nothing.
VERTEBRATE_MITOCHONDRIA = (
    ("MT", "Vertebrate_Mitochondrial"),
    ("M", "Vertebrate_Mitochondrial"),
)

# Keys are lowercase; look them up through mitochondria(), which folds case. config.yaml ships
# species capitalised ("Homo_sapiens"), so a case-sensitive lookup would miss every entry here.
MITOCHONDRIA_BY_SPECIES = {
    # NCBI table 5.
    "caenorhabditis_elegans": (("MtDNA", "Invertebrate_Mitochondrial"),),
    "drosophila_melanogaster": (("mitochondrion_genome", "Invertebrate_Mitochondrial"),),
    # NCBI table 3.
    "saccharomyces_cerevisiae": (("Mito", "Yeast_Mitochondrial"),),
    # NCBI table 13. Names its contig MT, so this is the one species the old hardcoded line reached
    # and mislabelled rather than missed.
    "ciona_intestinalis": (("MT", "Ascidian_Mitochondrial"),),
    # No mitochondrial contig in the Ensembl assembly, so nothing to declare.
    "ciona_savignyi": (),
}


def mitochondria(species):
    """The (contig, codon table) pairs to declare for a species, folding case on the name."""
    return MITOCHONDRIA_BY_SPECIES.get(species.lower(), VERTEBRATE_MITOCHONDRIA)


def snpeff_config_block(genome, species, assembly, release, gene_model=None, division="vertebrates"):
    """The text to append to snpEff.config for one genome, ending in a newline.

    `genome` is the name SnpEff will be asked to build, which is not always species.assembly - the
    isoform path builds a custom gene model under its own name. `gene_model` describes that when it
    applies, so the entry does not claim to be the plain Ensembl reference.

    `division` selects between Ensembl's main vertebrate site and Ensembl Bacteria, which differ in
    both the reference URL and the codon table.

    Note that SnpEff ships genome entries for bacteria already - 1086 Pseudomonas ones - but none of
    them carries a codonTable line, so every one of them translates with the standard table. The
    block written here is for a database Spritz builds itself, and does declare it.
    """
    description = f"{species} {assembly}"
    if gene_model:
        description = f"{description} with {gene_model}"

    reference = (ENSEMBL_BACTERIA_REFERENCE if division == "bacteria" else ENSEMBL_REFERENCE)

    lines = [
        "",
        f"# {genome}",
        f"{genome}.genome : {description}",
        f"{genome}.reference : {reference.format(release=release)}",
    ]
    if division == "bacteria":
        lines.append(f"\t{genome}.codonTable : {bacterial_codon_table(species)}")
    else:
        for contig, table in mitochondria(species):
            lines.append(f"\t{genome}.{contig}.codonTable : {table}")
    return "\n".join(lines) + "\n"
