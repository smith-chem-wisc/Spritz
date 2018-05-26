using System.Collections.Generic;

namespace WorkflowLayer
{
    public class SampleSpecificProteinDBParameters
        : ISpritzParameters
    {
        public string SpritzDirectory { get; set; }
        public string AnalysisDirectory { get; set; }
        public string Reference { get; set; }
        public int Threads { get; set; } = 1;
        public List<string[]> Fastqs { get; set; }
        public bool StrandSpecific { get; set; }
        public bool InferStrandSpecificity { get; set; }
        public bool OverwriteStarAlignment { get; set; }
        public string GenomeStarIndexDirectory { get; set; }
        public string GenomeFasta { get; set; }
        public string ProteinFasta { get; set; }
        public string ReferenceGeneModelGtfOrGff { get; set; }
        public string NewGeneModelGtfOrGff { get; set; }
        public string UniProtXmlPath { get; set; }
        public string EnsemblKnownSitesPath { get; set; }
        public int MinPeptideLength { get; set; } = 7;
        public bool UseReadSubset { get; set; } = false;
        public int ReadSubset { get; set; } = 300000;
        public bool DoFusionAnalysis { get; set; }
        public bool DoTranscriptIsoformAnalysis { get; set; }
    }
}