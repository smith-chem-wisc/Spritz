import os
import posixpath  # snakemake paths are always forward-slash; see the note on all_output

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
# Build output goes to a directory this workflow chooses, rather than to the SDK's default
# bin/<Platform>/<Configuration>/<TargetFramework>/ layout. Naming that default here meant the
# workflow silently depended on the platform, the configuration and the target framework all
# staying put, and would break as soon as any of them changed.
#
# Defined once, in two forms, because build_transfer_mods runs `dotnet build -o` from inside the
# project directory while the rest of the workflow refers to the assembly from the workflow
# directory. The rule passes SPRITZ_MODS_OUT_SUBDIR to -o via params, so the path the SDK writes to
# and the path declared as the rule's output cannot fall out of step.
SPRITZ_MODS_PROJECT_DIR = "../SpritzModifications"
SPRITZ_MODS_OUT_SUBDIR = "bin/spritz"
SPRITZ_MODS_OUT_DIR = f"{SPRITZ_MODS_PROJECT_DIR}/{SPRITZ_MODS_OUT_SUBDIR}"
TRANSFER_MOD_DLL="../SpritzModifications.dll" if PREBUILT_SPRITZ_MODS else f"{SPRITZ_MODS_OUT_DIR}/SpritzModifications.dll"

# .NET on Linux delegates TLS validation to OpenSSL, which resolves its CA store from its own
# compiled-in OPENSSLDIR and from SSL_CERT_FILE / SSL_CERT_DIR rather than from the active conda
# environment. Which store gets consulted therefore depends on which libssl happens to be loaded,
# and that varies between environments. Naming the bundle explicitly removes the ambiguity.
#
# Motivated by the report in issue #240, where the HTTPS download in setup_transfer_mods failed with
#
#   System.Security.Authentication.AuthenticationException: The remote certificate is invalid
#   because of errors in the certificate chain: UntrustedRoot
#     at UsefulProteomicsDatabases.Loaders.DownloadPsiMod
#
# while, in the same run on the same machine, git cloned over HTTPS from github.com and `dotnet
# restore` reached nuget.org over TLS - both successfully, in other conda environments. So the server
# certificate and the machine's trust generally were fine; only .NET inside this environment failed.
# The precise mechanism was never established and the report has not recurred, so treat this as
# removing a known source of variability rather than as a proven fix.
#
# Guarded on the file existing so that an environment already resolving trust correctly - notably the
# container, which works today - behaves exactly as before. Written as an `if` rather than
# `[ -f x ] && export ...` purely for legibility; both are safe under the bash strict mode snakemake
# uses, since a failing test that is not the final command of an && list does not trip errexit.
CONDA_CA_BUNDLE_EXPORT = (
    'if [ -f "$CONDA_PREFIX/ssl/cacert.pem" ]; '
    'then export SSL_CERT_FILE="$CONDA_PREFIX/ssl/cacert.pem"; fi; '
)

def all_output(wildcards):
    '''Gets the final output files depending on the configuration'''
    outputs = ["prose.txt"]
    outputs.append(posixpath.join("variants/", f"done{REF}.{ENSEMBL_VERSION}.txt")) # reference
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
        # posixpath, not os.path: on native Windows os.path.join inserts a backslash, so this
        # asked for a path no rule declares and the DAG failed before any job ran (issue #243).
        [posixpath.join("{dir}", file) for file in outputs],
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
