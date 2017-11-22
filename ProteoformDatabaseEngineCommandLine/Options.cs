using CommandLine;
using System;
using System.IO;
using System.Reflection;

namespace ProteoformDatabaseEngineCommandLine
{
    class Options
    {
        [Option('c', "command", Required = true, HelpText = "command: (1) setup, (2) run, (3) starFusionTest")]
        public string Command { get; set; }

        [Option('b', "binDirectory", Required = false, HelpText = "bin directory for PDE")]
        public string BinDirectory { get; set; } = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);

        [Option('a', "analysisDirectory", Required = false, HelpText = "Target directory for downloads and analysis")]
        public string AnalysisDirectory { get; set; } = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);

        [Option('r', "StarFusionReference", Required = false, HelpText = "Human reference for STAR fusion (GRCh37 or GRCh38)", DefaultValue = "GRCh38")]
        public string Reference { get; set; }

        [Option('t', "threads", Required = false, HelpText = "Number of threads to use")]
        public int Threads { get; set; } = Environment.ProcessorCount;

        [Option('s', "sraAccession", Required = false, HelpText = "SRR or SRX accession for download of fastq files for analysis")]
        public string SraAccession { get; set; }

        [Option("fq1", Required = false, HelpText = "FASTQ file for single-end or for pair1")]
        public string Fastq1 { get; set; }

        [Option("fq2", Required = false, HelpText = "FASTQ file pair2")]
        public string Fastq2 { get; set; }

        [Option('v', "overwriteStarAlignments", Required = false, HelpText = "overwrite STAR alignments if they already exist)")]
        public bool OverwriteStarAlignments { get; set; }

        [Option('s', "strandSpecific", Required = false, HelpText = "strandedness of the protocol", DefaultValue = false)]
        public bool StrandSpecific { get; set; }

        [Option('i', "inferStrandedness", Required = false, HelpText = "infer the strandedness with a sample of the fastq files", DefaultValue = true)]
        public bool InferStrandSpecificity { get; set; }

        [Option('d', "genomeDir", Required = true, HelpText = "STAR genome directory")]
        public string GenomeStarIndexDirectory { get; set; }

        [Option('f', "genomeFasta", Required = true, HelpText = "genomeFasta")]
        public string GenomeFasta { get; set; }

        [Option('g', "geneModelGtfOrGff", Required = true, HelpText = "geneModelGtfOrGff")]
        public string GeneModelGtfOrGff { get; set; }
    }
}
