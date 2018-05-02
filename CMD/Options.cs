using CommandLine;
using System;
using System.IO;
using System.Reflection;

namespace CMD
{
    internal class Options
    {
        [Option('c', "command", Required = true, HelpText = "Command: (1) setup, (2) run, (3) vcf2protein, (4) starFusionTest, (5) lncRNADiscovery, (6) quantify")]
        public string Command { get; set; }

        [Option('b', "binDirectory", Required = false, HelpText = "Bin directory for Proteoform Database Engine")]
        public string BinDirectory { get; set; } = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);

        [Option('a', "analysisDirectory", Required = false, HelpText = "Target directory for downloads and analysis")]
        public string AnalysisDirectory { get; set; } = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);

        [Option("fq1", Required = false, HelpText = "FASTQ file for single-end or for pair1 (comma-separated for multiple files)")]
        public string Fastq1 { get; set; }

        [Option("fq2", Required = false, HelpText = "FASTQ file pair2 (comma-separated for multiple files)")]
        public string Fastq2 { get; set; }

        [Option('s', "sraAccession", Required = false, HelpText = "SRR or SRX accession for download of fastq files for analysis (comma-separated for multiple SRAs)")]
        public string SraAccession { get; set; }

        [Option('t', "threads", Required = false, HelpText = "Number of threads to use")]
        public int Threads { get; set; } = Environment.ProcessorCount;

        [Option('d', "genomeDir", Required = false, HelpText = "STAR genome directory (default = genomeFastq without extension)")]
        public string GenomeStarIndexDirectory { get; set; }

        [Option('f', "genomeFasta", Required = false, HelpText = "Genome fasta (default = downloaded from Ensembl based on reference)")]
        public string GenomeFasta { get; set; }

        [Option('g', "geneModelGtfOrGff", Required = false, HelpText = "Gene model, either GTF, GFF2, or GFF3 (default = downloaded from Ensembl based on reference)")]
        public string GeneModelGtfOrGff { get; set; }

        [Option('v', "dbsnpVcfReference", Required = false, HelpText = "VCF reference file from dbSNP")]
        public string ReferenceVcf { get; set; }

        [Option('r', "StarFusionReference", Required = false, HelpText = "Human reference for STAR fusion (GRCh37 or GRCh38)", Default = "GRCh38")]
        public string Reference { get; set; }

        [Option('x', "UniProtProteinXml", Required = false, HelpText = "Protein XML UniProt Database for Homo sapiens")]
        public string UniProtXml { get; set; }

        [Option("overwriteStarAlignments", Required = false, HelpText = "Overwrite STAR alignments if they already exist", Default = false)]
        public bool OverwriteStarAlignments { get; set; }

        [Option("strandSpecific", Required = false, HelpText = "Stranded library preparation protocol", Default = false)]
        public bool StrandSpecific { get; set; }

        [Option("inferStrandedness", Required = false, HelpText = "Infer the strandedness with a sample of the FASTQ files", Default = false)]
        public bool InferStrandSpecificity { get; set; }
    }
}