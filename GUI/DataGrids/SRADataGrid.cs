namespace SpritzGUI
{
    internal class SRADataGrid
    {
        public SRADataGrid(string name)
        {
            Name = name;
        }

        public string Name { get; set; }
        public string State { get; set; }
        public int Completion { get; set; }
    }
}