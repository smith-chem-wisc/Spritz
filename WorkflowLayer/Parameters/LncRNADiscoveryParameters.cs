using System.Collections.Generic;

namespace WorkflowLayer
{
    public class LncRNADiscoveryParameters
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
        public string GeneModelGtfOrGff { get; set; }
        public bool DoOutputQuantificationBam { get; set; }
        public bool UseReadSubset { get; set; }
        public int ReadSubset { get; set; } = 300000;
    }
}