using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class TranscriptQuantificationParameters
        : ISpritzParameters
    {
        public TranscriptQuantificationParameters(string spritzDirectory, string referenceFastaPath, int threads, string geneModelPath, RSEMAlignerOption aligner,
            Strandedness strandedness, string[] fastq, bool doOutputBam)
        {
            SpritzDirectory = spritzDirectory;
            ReferenceFastaPath = referenceFastaPath;
            Threads = threads;
            GeneModelPath = geneModelPath;
            Aligner = aligner;
            Strandedness = strandedness;
            Fastq = fastq;
            DoOutputQuantificationBam = doOutputBam;
        }

        public string SpritzDirectory { get; }
        public string ReferenceFastaPath { get; }
        public int Threads { get; }
        public string GeneModelPath { get; }
        public RSEMAlignerOption Aligner { get; }
        public Strandedness Strandedness { get; }
        public bool DoOutputQuantificationBam { get; }
        public string[] Fastq { get; }
        public int ReadSubset { get; }
    }
}