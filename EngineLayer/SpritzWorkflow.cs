using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WorkflowLayer
{
    public enum MyWorkflow
    {
        Fastaq2Proteins,
        LncRnaDiscovery
    }
    public class SpritzWorkflow
    {
        protected SpritzWorkflow(MyWorkflow workflowType)
        {
            WorkflowType = workflowType;
        }

        public MyWorkflow WorkflowType { get; set; }

    }
}
