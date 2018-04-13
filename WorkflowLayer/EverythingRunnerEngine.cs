using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace WorkflowLayer
{
    public class EverythingRunnerEngine
    {
        private readonly List<Tuple<string, SpritzFlow>> taskList;
        private string outputFolder;
        private List<string> currentRnaSeqFilenameList;
        private List<string> currentGenomeFilenameList;
        private List<string> currentGeneSetFilenameList;

        public EverythingRunnerEngine(List<Tuple<string, SpritzFlow>> taskList, List<string> startingGenomeFilenameList, List<string> startingGeneSetFilenameList, List<string> startingRnaSeqFilenameList, string outputFolder)
        {
            this.taskList = taskList;
            this.outputFolder = outputFolder;

            currentGenomeFilenameList = startingGenomeFilenameList;
            currentGeneSetFilenameList = startingGeneSetFilenameList;
            currentRnaSeqFilenameList = startingRnaSeqFilenameList;
        }

        public static event EventHandler<StringListEventArgs> NewRnaSeqFastqHandler;

        public void Run()
        {
            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            //Used for test
            var writtenFile = Path.Combine(outputFolder, "test.txt");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("What a day!");
            }

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];

                var outputFolderForThisTask = Path.Combine(outputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                    Directory.CreateDirectory(outputFolderForThisTask);

                ok.Item2.RunTask(outputFolderForThisTask, currentGenomeFilenameList, currentGeneSetFilenameList, currentRnaSeqFilenameList, ok.Item1);
            }
        }

        private void NewRnaSeqFastq(List<string> newRnaSeqFastq)
        {
            NewRnaSeqFastqHandler?.Invoke(this, new StringListEventArgs(newRnaSeqFastq));
        }
    }
}