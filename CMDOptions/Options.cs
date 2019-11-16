using CommandLine;
using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using WorkflowLayer;

namespace CMD
{
    public class Options
    {
        public const char CommandOptionShort = 'c';
        public const string CommandOptionLong = "command";

        public string Command { get; set; }
        public string SpritzDirectory { get; set; } = Path.GetDirectoryName(Assembly.GetEntryAssembly().Location);
        public string AnalysisDirectory { get; set; } = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "output");
        public string Fastq1 { get; set; }
        public string Fastq2 { get; set; }
        public string ExperimentType { get; set; }
        public string SraAccession { get; set; }
        public int Threads { get; set; } = Environment.ProcessorCount == 1 ? 1 : Environment.ProcessorCount - 1;
        public string GenomeStarIndexDirectory { get; set; }
        public string GenomeFasta { get; set; }
        public string GeneModelGtfOrGff { get; set; }
        public string NewGeneModelGtfOrGff { get; set; }
        public string ReferenceVcf { get; set; }
        public string Reference { get; set; }
        public string UniProtXml { get; set; }
        public bool OverwriteStarAlignments { get; set; }
        public bool StrandSpecific { get; set; }
        public bool InferStrandSpecificity { get; set; }
        public bool SkipVariantAnalysis { get; set; }
        public bool DoTranscriptIsoformAnalysis { get; set; }
        public bool DoFusionAnalysis { get; set; }
        public string IndelFinder { get; set; }
        public int VariantCallingWorkers { get; set; } = 1;
        public string ProteinFastaPath { get; set; }

        // new for snakemake
        public string Release { get; set; }
        public string Species { get; set; }
        public string SnpEff { get; set; }
        public string Organism { get; set; }
    }
}