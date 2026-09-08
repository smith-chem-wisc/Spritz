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
        "Invertebrate_Mitochondrial",
        "Vertebrate_Mitochondrial",
        "Yeast_Mitochondrial",
    }
)

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


def snpeff_config_block(genome, species, assembly, release, gene_model=None):
    """The text to append to snpEff.config for one genome, ending in a newline.

    `genome` is the name SnpEff will be asked to build, which is not always species.assembly - the
    isoform path builds a custom gene model under its own name. `gene_model` describes that when it
    applies, so the entry does not claim to be the plain Ensembl reference.
    """
    description = f"{species} {assembly}"
    if gene_model:
        description = f"{description} with {gene_model}"

    lines = [
        "",
        f"# {genome}",
        f"{genome}.genome : {description}",
        f"{genome}.reference : https://ftp.ensembl.org/pub/release-{release}/",
    ]
    for contig, table in mitochondria(species):
        lines.append(f"\t{genome}.{contig}.codonTable : {table}")
    return "\n".join(lines) + "\n"
