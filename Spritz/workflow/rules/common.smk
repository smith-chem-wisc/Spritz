import os

# Variables used by many of the rules
SPECIES = config["species"]
GENOME_VERSION = config["genome"]
ENSEMBL_VERSION = config["release"]
GENEMODEL_VERSION = f"{GENOME_VERSION}.{ENSEMBL_VERSION}"
REF = f"{SPECIES}.{GENOME_VERSION}"
GENOME_FA = f"../resources/ensembl/{REF}.dna.primary_assembly.fa"
KARYOTYPIC_GENOME_PREFIX = f"../resources/ensembl/{REF}.dna.primary_assembly.karyotypic"
KARYOTYPIC_GENOME_FA = f"{KARYOTYPIC_GENOME_PREFIX}.fa"
ENSEMBL_GFF = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}.gff3"
TEST_GENOME_FA = f"../resources/ensembl/202122.fa"
TEST_ENSEMBL_GFF = f"../resources/ensembl/202122.gff3"
FA=GENOME_FA # for analysis; can also be TEST_GENOME_FA
GFF3=ENSEMBL_GFF # for analysis; can also be TEST_ENSEMBL_GFF
REFSTAR_PREFIX = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/RsemStarReference"
REFSTAR_FOLDER = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}RsemStar/"
REF_PREFIX = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/RsemReference"
REF_FOLDER = f"../resources/ensembl/{SPECIES}.{GENEMODEL_VERSION}Rsem/"
UNIPROTXML=f"../resources/uniprot/{config['species']}.protein.xml.gz" #"../resources/Homo_sapiens_202022.xml.gz"
UNIPROTFASTA=f"../resources/uniprot/{config['species']}.protein.fasta" #"../resources/Homo_sapiens_202022.xml.gz"
PREBUILT_SPRITZ_MODS = "prebuilt_spritz_mods" in config and config["prebuilt_spritz_mods"]
TRANSFER_MOD_DLL="../SpritzModifications.dll" if PREBUILT_SPRITZ_MODS else "../SpritzModifications/bin/x64/Release/net5.0/SpritzModifications.dll"

def all_output(wildcards):
    '''Gets the final output files depending on the configuration'''
    outputs = ["prose.txt"]
    else:
        outputs.append(os.path.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")) # reference
        if "variant" in config["analyses"]:
            outputs.append("final/combined.spritz.snpeff.protein.withmods.xml.gz") # variants
        if "isoform" in config["analyses"]:
            outputs.append("final/combined.spritz.isoform.protein.withmods.xml.gz") # isoforms
        if "quant" in config["analyses"]:
            outputs.extend([
                "final/transcript_reference_quant.tpms.csv",
                "final/gene_reference_quant.tpms.csv"]) # reference quant with stringtie
        if "variant" in config["analyses"] and "isoform" in config["analyses"]:
            outputs.append("final/combined.spritz.isoformvariants.protein.withmods.xml.gz") # isoform variants
        if "isoform" in config["analyses"] and "quant" in config["analyses"]:
            outputs.extend([
                "final/transcript_custom_quant.tpms.csv",
                "final/gene_custom_quant.tpms.csv"]) # isoform quant with stringtie
    expanded_outputs = expand(
        [os.path.join("{dir}", file) for file in outputs],
        dir=config["analysis_directory"])
    return expanded_outputs

def setup_output(wildcards):
    '''Gets the output needed for setting up spritz'''
    setup_outputs = [
        f"../resources/ChromosomeMappings/{GENOME_VERSION}_UCSC2ensembl.txt",
        f"../resources/ensembl/{SPECIES}.ensembl.vcf",
        TRANSFER_MOD_DLL,
        UNIPROTFASTA,
        UNIPROTFASTA,
        FA,
        "../resources/ptmlist.txt",
        "../resources/PSI-MOD.obo.xml",
        "../resources/SnpEff/snpEff.jar",
        [] if not check('sra') else expand(
            [
                "{dir}/{sra}_1.fastq",
                "{dir}/{sra}_2.fastq",
            ],
            dir=config['analysis_directory'],
            sra=config['sra']),
        [] if not check('sra_se') else expand(
            [
                "{dir}/{sra_se}.fastq",
            ],
            dir=config['analysis_directory'],
            sra_se=config['sra_se'])
    ]
    return setup_outputs

def check(field):
    '''Checks whether or not a field is contained in the configuration'''
    return field in config and config[field] is not None and len(config[field]) > 0
