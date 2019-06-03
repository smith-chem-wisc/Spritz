using CMD;
using Nett;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
using YamlDotNet.Core.Events;
using YamlDotNet.RepresentationModel;
using YamlDotNet.Serialization;

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

        public string Arguments { get; set; }
        public string StdErr { get; set; }
        public static string SpritzDirectory { get; set; } = Environment.CurrentDirectory;

        public void Run()
        {
            var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", System.Globalization.CultureInfo.InvariantCulture);

            outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];
                WriteConfig(ok.Item2);

                Process proc = new Process();
                proc.StartInfo.FileName = "Powershell.exe";
                proc.StartInfo.Arguments = "docker pull rinaibrhm/spritz ; docker run --rm -t -i -v \"\"\"" + ok.Item2.AnalysisDirectory + ":/app/data\"\"\" rinaibrhm/spritz";
                proc.StartInfo.CreateNoWindow = true;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardError = true;
                proc.Start();
                proc.WaitForExit();
                StdErr = proc.StandardError.ReadToEnd();
            }
        }

        public IEnumerable<string> GenerateCommandsDry()
        {
            for (int i = 0; i < taskList.Count; i++)
            {
                var options = taskList[i].Item2;
                yield return "docker pull rinaibrhm/spritz ; docker run --rm -t -i -v \"\"\"" + options.AnalysisDirectory + ":/app/data\"\"\" rinaibrhm/spritz";
            }
        }

        private string[] GenerateArguments(Options options)
        {
            List<string> commands = new List<string> { "-c", options.Command };
            if (options.SpritzDirectory != null && options.SpritzDirectory != "")
            {
                commands.AddRange(new[] { "-b", AddQuotes(options.SpritzDirectory) });
            }
            if (options.AnalysisDirectory != null && options.AnalysisDirectory != "")
            {
                commands.AddRange(new[] { "-a", AddQuotes(options.AnalysisDirectory)});
            }
            if (options.Fastq1 != null && options.Fastq1 != "")
            {
                commands.AddRange(new[] { "--fq1", AddQuotes(options.Fastq1) });
            }
            if (options.Fastq2 != null && options.Fastq2 != "")
            {
                commands.AddRange(new[] { "--fq2", AddQuotes(options.Fastq2) });
            }
            if (options.ExperimentType != null && options.ExperimentType != "")
            {
                commands.AddRange(new[] { "-e", options.ExperimentType });
            }
            if (options.SraAccession != null && options.SraAccession != "")
            {
                commands.AddRange(new[] { "-s", options.SraAccession });
            }
            if (options.Threads > 0 && options.Threads <= Environment.ProcessorCount)
            {
                commands.AddRange(new[] { "-t", options.Threads.ToString() });
            }
            if (options.GenomeStarIndexDirectory != null && options.GenomeStarIndexDirectory != "")
            {
                commands.AddRange(new[] { "-d", AddQuotes(options.GenomeStarIndexDirectory) });
            }
            if (options.GenomeFasta != null && options.GenomeFasta != "")
            {
                commands.AddRange(new[] { "-f", AddQuotes(options.GenomeFasta) });
            }
            if (options.GeneModelGtfOrGff != null && options.GeneModelGtfOrGff != "")
            {
                commands.AddRange(new[] { "-g", AddQuotes(options.GeneModelGtfOrGff) });
            }
            if (options.NewGeneModelGtfOrGff != null && options.NewGeneModelGtfOrGff != "")
            {
                commands.AddRange(new[] { "-h", AddQuotes(options.NewGeneModelGtfOrGff) });
            }
            if (options.ReferenceVcf != null && options.ReferenceVcf != "")
            {
                commands.AddRange(new[] { "-v", AddQuotes(options.ReferenceVcf) });
            }
            if (options.Reference != null && options.Reference != "")
            {
                commands.AddRange(new[] { "-r", options.Reference });
            }
            if (options.UniProtXml != null && options.UniProtXml != "")
            {
                commands.AddRange(new[] { "-x", AddQuotes(options.UniProtXml) });
            }
            if (options.UniProtXml != null && options.UniProtXml != "")
            {
                commands.AddRange(new[] { "--indelFinder", options.IndelFinder });
            }
            if (options.OverwriteStarAlignments)
            {
                commands.Add("--overwriteStarAlignments");
            }
            if (options.StrandSpecific)
            {
                commands.Add("--strandSpecific");
            }
            if (options.InferStrandSpecificity)
            {
                commands.Add("--inferStrandedness");
            }
            if (options.DoTranscriptIsoformAnalysis)
            {
                commands.Add("--doTranscriptIsoformAnalysis");
            }
            if (options.DoFusionAnalysis)
            {
                commands.Add("--doGeneFusionAnalysis");
            }
            if (options.SkipVariantAnalysis)
            {
                commands.Add("--skipVariantAnalysis");
            }
            if (options.VariantCallingWorkers > 0 && options.VariantCallingWorkers <= Environment.ProcessorCount)
            {
                commands.AddRange(new[] { "--variantCallingWorkers", options.VariantCallingWorkers.ToString() });
            }
            
            return commands.ToArray();
        }

        private void WriteConfig(Options options)
        {
            const string initialContent = "---\nversion: 1\n"; // needed to start writing yaml file

            var sr = new StringReader(initialContent);
            var stream = new YamlStream();
            stream.Load(sr);

            var rootMappingNode = (YamlMappingNode)stream.Documents[0].RootNode;

            var sras = options.SraAccession.Split(',');
            var fq1s = options.Fastq1.Split(',') ?? new string[0];
            var fq2s = options.Fastq2.Split(',')?? new string[0];

            // write user input sras
            var accession = new YamlSequenceNode();
            accession.Style = SequenceStyle.Flow;
            foreach (string sra in sras)
            {
                if (sra.Length > 0)
                {
                    accession.Add(sra);
                }
            }
            rootMappingNode.Add("sra", accession);

            // write user defined analysis directory (input and output folder)
            var analysisDirectory = new YamlSequenceNode();
            analysisDirectory.Style = SequenceStyle.Flow;
            analysisDirectory.Add(options.AnalysisDirectory);
            rootMappingNode.Add("analysisDirectory", analysisDirectory);

            // write user input fastqs
            var fq1 = new YamlSequenceNode();
            fq1.Style = SequenceStyle.Flow;
            foreach (string fastq in fq1s)
            {
                if (fastq.Length > 0)
                {
                    fq1.Add(fastq);
                }
            }
            rootMappingNode.Add("fq1", fq1);

            var fq2 = new YamlSequenceNode();
            fq2.Style = SequenceStyle.Flow;
            foreach (string fastq in fq2s)
            {
                if (fastq.Length > 0)
                {
                    fq2.Add(fastq);
                }
            }
            rootMappingNode.Add("fq2", fq2);

            Directory.SetCurrentDirectory(options.AnalysisDirectory); // switch to analysis directory for writing permissions
    
            // add config file to user defined analysis directory, will mount to docker container
            using (TextWriter writer = File.CreateText(Path.Combine(Directory.GetCurrentDirectory(), "config.yaml")))
            {
                stream.Save(writer, false);
            }

            // switch back to current directory
            Directory.SetCurrentDirectory(SpritzDirectory);
        }

        /// <summary>
        /// Needed for paths that have spaces in them to pass through properly
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        private static string AddQuotes(string path)
        {
            if (path.StartsWith("\"") && path.EndsWith("\""))
                return path;
            else
                return $"\"{path}\"";
        }
    }
}