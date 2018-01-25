using WorkflowLayer;

namespace SpritzGUI
{
    public class InRunTask : ForTreeView
    {
        #region Public Fields

        public readonly SpritzWorkflow workflow;

        #endregion Public Fields

        #region Public Constructors

        public InRunTask(string displayName, SpritzWorkflow workflow) : base(displayName, displayName)
        {
            this.workflow = workflow;
        }

        #endregion Public Constructors
    }
}