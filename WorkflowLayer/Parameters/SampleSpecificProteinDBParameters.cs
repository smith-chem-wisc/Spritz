using System.Collections.Generic;

namespace WorkflowLayer
{
    public class SampleSpecificProteinDBParameters
        : ISpritzParameters
    {
        public SampleSpecificProteinDBParameters(string spritzDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath,
            bool useReadSubset = false, int readSubset = 300000)
        {
            SpritzDirectory = spritzDirectory;
            AnalysisDirectory = analysisDirectory;
            Reference = reference;
            Threads = threads;
            Fastqs = fastqs;
            StrandSpecific = strandSpecific;
            InferStrandSpecificity = inferStrandSpecificity;
            OverwriteStarAlignment = overwriteStarAlignment;
            GenomeStarIndexDirectory = genomeStarIndexDirectory;
            GenomeFasta = genomeFasta;
            ProteinFasta = proteinFasta;
            GeneModelGtfOrGff = geneModelGtfOrGff;
            EnsemblKnownSitesPath = ensemblKnownSitesPath;
            UseReadSubset = useReadSubset;
            ReadSubset = readSubset;
        }

        public string SpritzDirectory { get; }
        public string AnalysisDirectory { get; }
        public string Reference { get; }
        public int Threads { get; }
        public List<string[]> Fastqs { get; }
        public bool StrandSpecific { get; }
        public bool InferStrandSpecificity { get; }
        public bool OverwriteStarAlignment { get; }
        public string GenomeStarIndexDirectory { get; }
        public string GenomeFasta { get; }
        public string ProteinFasta { get; }
        public string GeneModelGtfOrGff { get; }
        public string EnsemblKnownSitesPath { get; }
        public bool UseReadSubset { get; }
        public int ReadSubset { get; } = 300000;
    }
}