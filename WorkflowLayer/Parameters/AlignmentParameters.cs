using ToolWrapperLayer;
using System.Collections.Generic;

namespace WorkflowLayer
{
    public class AlignmentParameters
        : ISpritzParameters
    {
        public AlignmentParameters()
        {
            Reference = "GRCh38";
            Threads = 12;
            StrandSpecific = false;
            InferStrandSpecificity = false;
            OverwriteStarAlignment = false;
            UseReadSubset = false;
            ReadSubset = 300000;
        }

        public string SpritzDirectory { get; set; }
        public string AnalysisDirectory { get; set; }
        public string Reference { get; set; }
        public int Threads { get; set; } = 1;
        public List<string[]> Fastqs { get; set; }
        public ExperimentType ExperimentType { get; set; }
        public bool StrandSpecific { get; set; }
        public bool InferStrandSpecificity { get; set; }
        public bool OverwriteStarAlignment { get; set; }
        public string GenomeStarIndexDirectory { get; set; }
        public string ReorderedFastaPath { get; set; }
        public string GeneModelGtfOrGffPath { get; set; }
        public string EnsemblKnownSitesPath { get; set; }
        public List<string> FirstPassSpliceJunctions { get; set; }
        public string SecondPassGenomeDirectory { get; set; }
        public List<string> SortedBamFiles { get; set; }
        public List<string> DedupedBamFiles { get; set; }
        public List<string> ChimericSamFiles { get; set; }
        public List<string> ChimericJunctionFiles { get; set; }
        public bool UseReadSubset { get; set; }
        public int ReadSubset { get; set; } = 300000;
    }
}