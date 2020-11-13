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

        public bool Use { get; set; }
        public string FileName { get; set; }
        public bool IsPairedEnd { get; set; }
        public string MatePair { get; set; }
        public string Experiment { get; set; }
        public string FilePath { get; set; }
    }
}