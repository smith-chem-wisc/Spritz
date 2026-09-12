import sys

import spritz_config

# From the rule via SPRITZ_CONFIG. Reading the relative path read the packaged defaults, so the
# methods text described a Homo_sapiens quant run whatever was actually run.
config = spritz_config.load()

workflows = config["analyses"]
outf = sys.argv[1]

# SPRITZ
lines = [
    "Spritz is a program that generates protein databases annotated with sequence variations and PTMs: ",
    "- Cesnik, A. J.; et al. Spritz: A Proteogenomic Database Engine. J. Proteome Res. 2020, in press. https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.0c00407.",
    "",
]

def supplied(field):
    """The same test check() makes in common.smk, so citations match the rules that ran."""
    return spritz_config.check(config, field)

# A supplied VCF replaces variant *calling*, but not alignment: quant and isoform reconstruction
# still align reads, and all_output asks for them independently of the VCF. So the trimming and
# alignment citations depend on whether reads were supplied, and only the GATK calling citation
# depends on the VCF. Citing tools that did not run, and omitting tools that did, are both wrong.
used_vcf = supplied("vcf")
used_reads = any(supplied(f) for f in ("sra", "sra_se", "fq", "fq_se"))

# SRA DOWNLOADS
used_sras = supplied("sra")
if used_sras:
    lines.extend([
        "SRAs are downloaded using the SRA toolkit from NCBI: ",
        "- Leinonen, R.; et al. International Nucleotide Sequence Database Collaboration. The Sequence Read Archive. Nucleic Acids Res. 2011, 39 (Database issue), D19-21. https://doi.org/10.1093/nar/gkq1019.",
        ""
    ])

# TRIMMING AND ALIGNMENT
if used_reads:
    lines.extend([
        "Reads are trimmed and analyzed for quality scores using fastp: ",
        "- Chen, S.; et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics 2018, 34 (17), i884-i890. https://academic.oup.com/bioinformatics/article/34/17/i884/5093234.",
        "",
        "Reads are aligned using hisat2: ",
        "- Kim, D.; et al. Graph-Based Genome Alignment and Genotyping with HISAT2 and HISAT-Genotype. Nat. Biotechnol. 2019, 37 (8), 907-915. https://doi.org/10.1038/s41587-019-0201-4.",
        "",
        "Alignments are analyzed and combined using samtools: ",
        "- Li, H.; et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 2009, 25 (16), 2078-2079. https://academic.oup.com/bioinformatics/article/25/16/2078/204688.",
        ""
    ])

# VARIANT CALLING
if "variant" in workflows:
    if used_vcf:
        lines.extend([
            f"Variants were supplied as {', '.join(config['vcf'])}, called outside this workflow, "
            f"rather than called from reads"
            + (f", and merged across {len(config['vcf'])} per-sample VCFs. " if len(config['vcf']) > 1 else ". "),
            ""])
    else:
        lines.extend([
            "Alignments are prepared for variant calling and analyzed for variants using the Genome Analysis Toolkit (GATK): ",
            "- McKenna, A.; et al. The Genome Analysis Toolkit: A MapReduce Framework for Analyzing next-Generation DNA Sequencing Data. Genome Res. 2010, 20 (9), 1297-1303. https://doi.org/10.1101/gr.107524.110.",
            "- DePristo, M. A.; et al. A Framework for Variation Discovery and Genotyping Using Next-Generation DNA Sequencing Data. Nat. Genet. 2011, 43 (5), 491-498. https://doi.org/10.1038/ng.806.",
            "- Poplin, R.; et al. Scaling Accurate Genetic Variant Discovery to Tens of Thousands of Samples; preprint; Genomics, 2017. https://doi.org/10.1101/201178.",
            ""])
        # Which known-set the recalibration used is a methods detail, and for most species it is
        # not a published one. Saying so is the difference between a reproducible methods section
        # and one that implies a reference set that does not exist.
        if (config.get("known_sites") or "ensembl") == "bootstrap":
            lines.extend([
                "No published set of known variant sites exists for this species, so base quality "
                "score recalibration used a bootstrapped set: variants were called once on "
                "unrecalibrated alignments, hard-filtered to high-confidence SNPs (GATK RNA-seq "
                "filters, QUAL >= 30), and used as the known sites for recalibration before the "
                "reported calling pass. ",
                ""])
    lines.extend([
        "SnpEff is used for variant annotation and customized in Spritz to output a proteogenomic database: ",
        "- Cingolani, P.; et al. A Program for Annotating and Predicting the Effects of Single Nucleotide Polymorphisms, SnpEff: SNPs in the Genome of Drosophila Melanogaster Strain W1118; Iso-2; Iso-3. Fly (Austin) 2012, 6 (2), 80-92. https://doi.org/10.4161/fly.19695.",
        ""])

# ISOFORM ANALYSIS
if "isoform" in workflows:
    lines.extend([
        "Isoform analysis is performed using a pipeline from ProteomeGenerator: ",
        "- Cifani, P.; et al. ProteomeGenerator: A Framework for Comprehensive Proteomics Based on de Novo Transcriptome Assembly and High-Accuracy Peptide Mass Spectral Matching. J. Proteome Res. 2018, 17 (11), 3681-3692. https://doi.org/10.1021/acs.jproteome.8b00295.",
        ""])

transcripts = "reference" if not "isoform" in workflows else "reference and assembled"
if "quant" in workflows:
    lines.extend([
        f"Transcript quantification was performed for {transcripts} transcripts using StringTie2: ",
        "- Kovaka, S.; et al. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biol 2019, 20 (278), 1-13. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1",
        ""])

with open(outf, 'w') as file:
    file.writelines([f"{x}\n" for x in lines]) # note: all text above needs to be UTF-8 characters
