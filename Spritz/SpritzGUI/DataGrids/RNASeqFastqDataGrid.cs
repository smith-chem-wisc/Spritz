using System.IO;

namespace Spritz
{
    internal class RNASeqFastqDataGrid
    {
        public RNASeqFastqDataGrid(string filePath, bool isPairedEnd)
        {
            Use = true;
            FilePath = filePath;
            FileName = Path.GetFileNameWithoutExtension(filePath);
            IsPairedEnd = isPairedEnd;
            if (filePath.EndsWith("gz"))
                FileName = Path.GetFileNameWithoutExtension(FileName);
            if (FileName.EndsWith("_1"))
                MatePair = "1";
            if (FileName.EndsWith("_2"))
                MatePair = "2";
        }

        /// <summary>
        /// A VCF called elsewhere, listed here rather than in a box of its own because reads and
        /// variant calls are alternative entry points to the same pipeline, so seeing them in one
        /// list is what makes it obvious they cannot both be used.
        /// </summary>
        public static RNASeqFastqDataGrid Vcf(string filePath) =>
            new(filePath, isPairedEnd: false) { IsVcf = true };

        public bool Use { get; set; }
        public string FileName { get; set; }
        public bool IsPairedEnd { get; set; }
        public string MatePair { get; set; }
        public string Experiment { get; set; }
        public string FilePath { get; set; }
        public bool IsVcf { get; set; }

        /// <summary>Shown as a column, since the grid generates one per public property.</summary>
        public string Kind => IsVcf ? "VCF" : "FASTQ";

        /// <summary>
        /// The name the workflow is given. FASTQs go in as a prefix with the mate suffix removed,
        /// which is what the align rules expect; a VCF is named in full, extension included.
        /// </summary>
        public string WorkflowName => IsVcf
            ? Path.GetFileName(FilePath)
            : FileName.Substring(0, FileName.Length - 2);
    }
}