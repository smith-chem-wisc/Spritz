using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WorkflowLayer
{
    public class EverythingRunnerEngine
    {
        private readonly List<Tuple<string, SpritzWorkflow>> taskList;
        private string outputFolder;
        List<string> currentRnaSeqFilenameList;
        List<string> currentGenomeFilenameList;
        List<string> currentGeneSetFilenameList;

        public EverythingRunnerEngine(List<Tuple<string, SpritzWorkflow>> taskList, List<string> startingRnaSeqFilenameList, List<string> startingGenomeFilenameList, List<string> startingGeneSetFilenameList, string outputFolder)
        {
            this.taskList = taskList;
            this.outputFolder = outputFolder;

            currentRnaSeqFilenameList = startingRnaSeqFilenameList;
            currentGenomeFilenameList = startingGenomeFilenameList;
            currentGeneSetFilenameList = startingGeneSetFilenameList;
        }

        public void Run()
        {

        }
    }
}
