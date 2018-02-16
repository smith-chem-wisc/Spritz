using WorkflowLayer;

namespace SpritzGUI
{
    internal class PreRunTask
    {
        #region Public Fields

        public readonly SpritzFlow spritzWorkflow;

        #endregion Public Fields

        #region Public Constructors

        public PreRunTask(SpritzFlow theWorkflow)
        {
            this.spritzWorkflow = theWorkflow;
        }

        #endregion Public Constructors

        #region Public Properties

        public string DisplayName { get; set; }

        #endregion Public Properties
    }
}