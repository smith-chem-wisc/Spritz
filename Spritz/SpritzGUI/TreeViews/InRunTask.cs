using SpritzBackend;

namespace Spritz
{
    public class InRunTask : ForTreeView
    {
        public readonly SpritzOptions options;

        public InRunTask(string displayName, SpritzOptions options) : base(displayName, displayName)
        {
            this.options = options;
        }
    }
}