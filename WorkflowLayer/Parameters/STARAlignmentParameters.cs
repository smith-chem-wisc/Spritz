using System.Collections.Generic;
using System.IO;

namespace WorkflowLayer
{
    public class STARAlignmentParameters
        : ISpritzParameters
    {
        public STARAlignmentParameters(string spritzDirectory, string analysisDirectory, string reference, int threads,
            List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string reorderedFasta, string geneModelGtfOrGff,
            bool useReadSubset = false, int readSubset = 30000)
        {
            SpritzDirectory = spritzDirectory;
            AnalysisDirectory = analysisDirectory;
            Reference = reference;
            Threads = threads;
            Fastqs = fastqs;
            StrandSpecific = strandSpecific;
            InferStrandSpecificity = inferStrandSpecificity;
            OverWriteStarAlignment = overwriteStarAlignment;
            GenomeStarIndexDirectory = genomeStarIndexDirectory;
            ReorderedFasta = reorderedFasta;
            GeneModelGtfOrGff = geneModelGtfOrGff;
            UseReadSubset = useReadSubset;
            ReadSubset = readSubset;
        }

        public STARAlignmentParameters()
        {
            Reference = "GRCh38";
            Threads = 12;
            StrandSpecific = false;
            InferStrandSpecificity = false;
            OverWriteStarAlignment = false;
            UseReadSubset = false;
            ReadSubset = 300000;
        }

        public string SpritzDirectory { get; set; }
        public string AnalysisDirectory { get; set; }
        public string Reference { get; set; }
        public int Threads { get; set; } = 1;
        public List<string[]> Fastqs { get; set; }
        public bool StrandSpecific { get; set; }
        public bool InferStrandSpecificity { get; set; }
        public bool OverWriteStarAlignment { get; set; }
        public string GenomeStarIndexDirectory { get; set; }
        public string ReorderedFasta { get; set; }
        public string GeneModelGtfOrGff { get; set; }
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