using WorkflowLayer;

namespace SpritzGUI
{
    internal class PreRunTask
    {
        #region Public Fields

        public readonly SpritzWorkflow spritzWorkflow;

        #endregion Public Fields

        #region Public Constructors

        public PreRunTask(SpritzWorkflow theWorkflow)
        {
            this.spritzWorkflow = theWorkflow;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; set; }

        #endregion Public Properties
    }
}