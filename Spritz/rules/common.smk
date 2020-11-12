import os
import numpy as np

SPECIES = config["species"]
GENOME_VERSION = config["genome"]
ENSEMBL_VERSION = config["release"]
GENEMODEL_VERSION = f"{GENOME_VERSION}.{ENSEMBL_VERSION}"
GENOME_FA = f"data/ensembl/{SPECIES}.{GENOME_VERSION}.dna.primary_assembly.fa"
ENSEMBL_GFF = f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}.gff3"
TEST_GENOME_FA = f"data/ensembl/202122.fa"
TEST_ENSEMBL_GFF = f"data/ensembl/202122.gff3"
FA=GENOME_FA # for analysis; can also be TEST_GENOME_FA
GFF3=ENSEMBL_GFF # for analysis; can also be TEST_ENSEMBL_GFF
REFSTAR_PREFIX = f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/"
REF_PREFIX = f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/RsemReference"
REF_FOLDER = f"data/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/"
REF = f"{SPECIES}.{GENOME_VERSION}"
UNIPROTXML=f"data/uniprot/{config['species']}.protein.xml.gz" #"data/Homo_sapiens_202022.xml.gz"
UNIPROTFASTA=f"data/uniprot/{config['species']}.protein.fasta" #"data/Homo_sapiens_202022.xml.gz"

def all_output(wildcards):
    outputs = []
    if not "spritzversion" in config:
        outputs = ["please_update_spritz.txt"]
    elif len(config["analyses"]) == 0:
        outputs = ["prose.txt",
            os.path.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")] # reference
    elif "variant" in config["analyses"] and len(config["analyses"]) == 1:
        outputs = ["prose.txt",
            "final/combined.spritz.snpeff.protein.withmods.xml.gz", # variants
            os.path.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")] # reference
    elif "isoform" in config["analyses"] and len(config["analyses"]) == 1:
        outputs = ["prose.txt",
            "final/combined.spritz.isoform.protein.withmods.xml.gz"] # isoforms
    elif "variant" in config["analyses"] and "isoform" in config["analyses"]:
        outputs = ["prose.txt",
            "final/combined.spritz.snpeff.protein.withmods.xml.gz", # variants
            "final/combined.spritz.isoformvariants.protein.withmods.xml.gz", # isoform variants
            "final/combined.spritz.isoform.protein.withmods.xml.gz", # isoforms
            os.path.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")] # reference
    expanded_outputs = expand(
        [os.path.join("{dir}", file) for file in outputs],
        dir=config["analysisDirectory"])
    return expanded_outputs

def check(field):
    return field in config and config[field] is not None and len(config[field]) > 0
