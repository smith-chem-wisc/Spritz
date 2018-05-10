using System.Collections.Generic;

namespace WorkflowLayer
{
    public class SampleSpecificProteinDBParameters
        : ISpritzParameters
    {
        public SampleSpecificProteinDBParameters(string spritzDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string referenceGeneModelGtfOrGff, string ensemblKnownSitesPath, string uniprotXmlPath,
            string newGeneModelGtfOrGff = null, 
            double readDepthPerMilliionVariantCutoff = 0.05, int minPeptideLength = 7, bool useReadSubset = false, int readSubset = 300000)
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
            ReferenceGeneModelGtfOrGff = referenceGeneModelGtfOrGff;
            NewGeneModelGtfOrGff = newGeneModelGtfOrGff;
            UniProtXmlPath = uniprotXmlPath;
            EnsemblKnownSitesPath = ensemblKnownSitesPath;
            ReadDepthPerMilliionVariantCutoff = readDepthPerMilliionVariantCutoff;
            MinPeptideLength = minPeptideLength;
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
        public string ReferenceGeneModelGtfOrGff { get; }
        public string NewGeneModelGtfOrGff { get; }
        public string UniProtXmlPath { get; }
        public string EnsemblKnownSitesPath { get; }
        public double ReadDepthPerMilliionVariantCutoff { get; }
        public int MinPeptideLength { get; }
        public bool UseReadSubset { get; }
        public int ReadSubset { get; } = 300000;
    }
}