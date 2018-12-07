using CMD;

namespace SpritzGUI
{
    public class InRunTask : ForTreeView
    {
        public readonly Options options;

        public InRunTask(string displayName, Options options) : base(displayName, displayName)
        {
            this.options = options;
        }
    }
}