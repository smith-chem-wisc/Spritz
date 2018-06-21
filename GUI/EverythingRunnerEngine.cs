using CMD;
using System;
using System.Collections.Generic;
using System.IO;

namespace SpritzGUI
{
    public class EverythingRunnerEngine
    {
        private readonly List<Tuple<string, Options>> taskList;
        private string outputFolder;

        public EverythingRunnerEngine(List<Tuple<string, Options>> taskList, string outputFolder)
        {
            this.taskList = taskList;
            this.outputFolder = outputFolder;
        }

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

                //Put it into a function
                Spritz.Main(new string[] { "CMD.exe", "-c", ok.Item2.Command,
                "b", ok.Item2.SpritzDirectory,
                "a", ok.Item2.AnalysisDirectory,
                "fq1", ok.Item2.Fastq1,
                "fq2", ok.Item2.Fastq2,
                "s", ok.Item2.SraAccession,
                "t", ok.Item2.Threads.ToString(),
                "d", ok.Item2.GenomeStarIndexDirectory,
                "f", ok.Item2.GenomeFasta,
                "g", ok.Item2.GeneModelGtfOrGff,
                "v", ok.Item2.ReferenceVcf,
                "r", ok.Item2.Reference,
                "x", ok.Item2.UniProtXml,
                "overwriteStarAlignments", ok.Item2.OverwriteStarAlignments.ToString(),
                "strandSpecific", ok.Item2.StrandSpecific.ToString(),
                "inferStrandedness", ok.Item2.InferStrandSpecificity.ToString()});
            }
        }
    }
}