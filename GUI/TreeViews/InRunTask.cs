using WorkflowLayer;

namespace SpritzGUI
{
    public class InRunTask : ForTreeView
    {
        #region Public Fields

        public readonly SpritzFlow workflow;

        #endregion Public Fields

        #region Public Constructors

        public InRunTask(string displayName, SpritzFlow workflow) : base(displayName, displayName)
        {
            this.workflow = workflow;
        }

        #endregion Public Constructors
    }
}