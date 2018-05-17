using System;
using System.Diagnostics;
using System.IO;

namespace WorkflowLayer
{
    public enum MyWorkflow
    {
        SampleSpecificProteinDB,
        LncRnaDiscovery,
        TranscriptQuantification,
        STARAlignment,
    }

    public abstract class SpritzFlow
    {
        protected SpritzFlow(MyWorkflow workflowType)
        {
            WorkflowType = workflowType;
        }

        public MyWorkflow WorkflowType { get; set; }

        public ISpritzParameters SpritzParameters { get; set; }

        public void RunTask(string outputFolder, ISpritzParameters parameters, string displayName)
        {
            try
            {
                var stopWatch = new Stopwatch();
                stopWatch.Start();
                RunSpecific(parameters);
                stopWatch.Stop();
                var resultsFileName = Path.Combine(outputFolder, "results.txt");
                using (StreamWriter file = new StreamWriter(resultsFileName))
                {
                    file.WriteLine("Spritz: version ");
                    file.Write(stopWatch.Elapsed.ToString());
                }
            }
            catch (Exception)
            {
                throw;
            }
        }

        protected abstract void RunSpecific(ISpritzParameters parameters);
    }
}