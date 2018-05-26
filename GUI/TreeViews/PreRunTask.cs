using WorkflowLayer;
using CMD;

namespace SpritzGUI
{
    internal class PreRunTask
    {
        #region Public Fields

        public readonly Options options;

        #endregion Public Fields

        #region Public Constructors

        public PreRunTask(Options options)
        {
            this.options = options;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; set; }

        #endregion Public Properties
    }
}