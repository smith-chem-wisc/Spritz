namespace SpritzGUI
{
    internal class GeneSetDataGrid
    {
        public GeneSetDataGrid(string filePath)
        {
            Use = true;
            FilePath = filePath;
        }

        public bool Use { get; set; }
        public bool InProgress { get; private set; }
        public string FilePath { get; set; }

        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }
    }
}