using CMD;

namespace SpritzGUI
{
    public class InRunTask : ForTreeView
    {
        #region Public Fields

        public readonly Options options;

        #endregion Public Fields

        #region Public Constructors

        public InRunTask(string displayName, Options options) : base(displayName, displayName)
        {
            this.options = options;
        }

        #endregion Public Constructors
    }
}