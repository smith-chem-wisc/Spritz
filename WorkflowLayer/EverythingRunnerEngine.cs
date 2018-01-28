using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;
using System.IO;


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

        public static event EventHandler<StringListEventArgs> NewRnaSeqFastqHandler;

        public void Run()
        {
            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            //Used for test
            //var writtenFile = Path.Combine(outputFolder, "test.txt");
            //using (StreamWriter output = new StreamWriter(writtenFile))
            //{
            //    output.WriteLine("What a day!");
            //}

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];

                var outputFolderForThisTask = Path.Combine(outputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);



            }
        }

        private void NewRnaSeqFastq(List<string> newRnaSeqFastq)
        {
            NewRnaSeqFastqHandler?.Invoke(this, new StringListEventArgs(newRnaSeqFastq));
        }
    }
}
