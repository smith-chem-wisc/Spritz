using System.IO;

namespace SpritzGUI
{
    internal class RNASeqFastqDataGrid
    {
        public RNASeqFastqDataGrid(string filePath)
        {
            Use = true;
            FileName = Path.GetFileNameWithoutExtension(filePath);
            FilePath = filePath;
        }

        public bool Use { get; set; }
        public string FileName { get; set; }
        public string Experiment { get; set; }
        public string MateRun { get; set; }
        public string FilePath { get; set; }

        public void SetConditionText(string text)
        {
            Experiment = text;
        }
        public void SetTechRepText(string text)
        {
            MateRun = text;
        }
    }
}