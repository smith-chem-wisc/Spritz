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

        public EverythingRunnerEngine(List<Tuple<string, SpritzFlow>> taskList, string outputFolder)
        {
            this.taskList = taskList;
            this.outputFolder = outputFolder;
        }

        public static event EventHandler<StringListEventArgs> NewRnaSeqFastqHandler;

        public void Run()
        {
            //var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            //outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

#if DEBUG
            var writtenFile = Path.Combine(outputFolder, "test.txt");

            if (!Directory.Exists(outputFolder))
                Directory.CreateDirectory(outputFolder);

            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine("What a wonderful day, what a wonderful Spritz!");
            }

#endif

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];

                var outputFolderForThisTask = Path.Combine(outputFolder, ok.Item1);

                if (!Directory.Exists(outputFolderForThisTask))
                {
                    Directory.CreateDirectory(outputFolderForThisTask);
                }                

                ok.Item2.RunTask(outputFolderForThisTask, ok.Item2.SpritzParameters ,ok.Item1);
            }
        }

    }
}