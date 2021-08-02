namespace Spritz
{
    internal class SRADataGrid
    {
        public SRADataGrid(string name, bool isPairedEnd)
        {
            Name = name;
            IsPairedEnd = isPairedEnd;
        }

        public string Name { get; set; }
        public bool IsPairedEnd { get; set; }
    }
}