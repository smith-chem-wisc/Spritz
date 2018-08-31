using CMD;
using System;
using System.Collections.Generic;
using System.IO;
using Nett;

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

                var tomlFileName = Path.Combine(outputFolder, i.ToString() + "_Parameters.toml");
                Toml.WriteFile(ok.Item2, tomlFileName);

                //Put it into a function
                var commands = generateCommand(ok.Item2);
                Spritz.Main(commands);
            }
        }

        private string[] generateCommand(Options options)
        {
            List<string> commands = new List<string>();
            //commands.Add("CMD.exe");
            commands.Add("-c");
            commands.Add(options.Command);
            if (options.SpritzDirectory != "")
            {
                commands.Add("b");
                commands.Add(options.SpritzDirectory);
            }
            if (options.AnalysisDirectory != "")
            {
                commands.Add("a");
                commands.Add(options.AnalysisDirectory);
            }
            if (options.Fastq1 != "" && options.Fastq1!=null)
            {
                commands.Add("fq1");
                commands.Add(options.Fastq1);
            }
            if (options.Fastq2!="" && options.Fastq2!=null)
            {
                commands.Add("fq2");
                commands.Add(options.Fastq2);
            }
            if (options.SraAccession != "")
            {
                commands.Add("s");
                commands.Add(options.SraAccession);
            }
            if (options.Threads > 0)
            {
                commands.Add("t");
                commands.Add(options.Threads.ToString());
            }
            if (options.GenomeStarIndexDirectory != "")
            {
                commands.Add("d");
                commands.Add(options.GenomeStarIndexDirectory);
            }
            if (options.GenomeFasta != "")
            {
                commands.Add("f");
                commands.Add(options.GenomeFasta);
            }
            if (options.GeneModelGtfOrGff != "")
            {
                commands.Add("g");
                commands.Add(options.GeneModelGtfOrGff);
            }
            if (options.NewGeneModelGtfOrGff != "")
            {
                commands.Add("h");
                commands.Add(options.NewGeneModelGtfOrGff);
            }
            if (options.ReferenceVcf != "")
            {
                commands.Add("v");
                commands.Add(options.ReferenceVcf);
            }
            if (options.Reference != "")
            {
                commands.Add("r");
                commands.Add(options.Reference);
            }
            if (options.UniProtXml != "")
            {
                commands.Add("x");
                commands.Add(options.UniProtXml);
            }
            if (options.OverwriteStarAlignments == true)
            {
                commands.Add("overwriteStarAlignments");
                commands.Add("true");
            }
            if (options.StrandSpecific == true)
            {
                commands.Add("strandSpecific");
                commands.Add("true");
            }
            if (options.InferStrandSpecificity == true)
            {
                commands.Add("inferStrandedness");
                commands.Add("true");
            }
            if (options.DoTranscriptIsoformAnalysis == true)
            {
                commands.Add("doTranscriptIsoformAnalysis");
                commands.Add("true");
            }
            if (options.DoFusionAnalysis==true)
            {
                commands.Add("doGeneFusionAnalysis");
                commands.Add("true");
            }
            if (options.QuickSnpEffWithoutStats == true)
            {
                commands.Add("quickSnpEffWithoutStats");
                commands.Add("true");
            }
            return commands.ToArray();
        }
    }
}