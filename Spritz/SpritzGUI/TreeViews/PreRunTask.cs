using SpritzBackend;

namespace Spritz
{
    internal class PreRunTask
    {
        public readonly SpritzOptions options;

        public PreRunTask(SpritzOptions options)
        {
            this.options = options;
        }

        public string DisplayName { get; set; }
    }
}