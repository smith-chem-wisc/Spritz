using CMD;
using Nett;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using YamlDotNet.Core.Events;
using YamlDotNet.RepresentationModel;

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
        //public static string SpritzDirectory { get; set; } = Environment.CurrentDirectory;
        public static string AnalysisDirectory { get; set; }
        public static string ConfigDirectory { get; set; }

        public void Run()
        {
            //var startTimeForAllFilenames = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture);

            //outputFolder = outputFolder.Replace("$DATETIME", startTimeForAllFilenames);

            for (int i = 0; i < taskList.Count; i++)
            {
                var ok = taskList[i];
                WriteConfig(ok.Item2);

                Process proc = new Process();
                proc.StartInfo.FileName = "Powershell.exe";
                proc.StartInfo.Arguments = "docker pull rinaibrhm/spritz ; docker run --rm -t -i --name spritz -v \"\"\"" + ok.Item2.AnalysisDirectory + ":/app/" + AnalysisDirectory + "\"\"\" -v \"\"\"" + ConfigDirectory + ":/app/configs\"\"\" rinaibrhm/spritz";
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

        private YamlSequenceNode AddParam(string[] items, YamlSequenceNode node)
        {
            node.Style = SequenceStyle.Flow;
            foreach (string item in items)
            {
                if (item.Length > 0)
                {
                    node.Add(item);
                }
            }

            return node;
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
            AnalysisDirectory = options.AnalysisDirectory.Split('\\').ToList().Last();

            // write user input sras
            var accession = new YamlSequenceNode();
            rootMappingNode.Add("sra", AddParam(sras, accession));

            // write user defined analysis directory (input and output folder)
            var analysisDirectory = new YamlSequenceNode();
            analysisDirectory.Style = SequenceStyle.Flow;
            analysisDirectory.Add(AnalysisDirectory);
            rootMappingNode.Add("analysisDirectory", analysisDirectory);

            // write user input fastqs
            var fq1 = new YamlSequenceNode();
            rootMappingNode.Add("fq1", AddParam(fq1s, fq1));
            var fq2 = new YamlSequenceNode();
            rootMappingNode.Add("fq2", AddParam(fq2s, fq2));
            
            var pathToConfig = Path.Combine(Directory.GetCurrentDirectory(), "configs");

            // create configs folder to mount newly written config to container
            if (!Directory.Exists(pathToConfig))
            {
                Directory.CreateDirectory(pathToConfig);
            }
            ConfigDirectory = pathToConfig;


            using (TextWriter writer = File.CreateText(Path.Combine(pathToConfig, "config.yaml")))
            {
                stream.Save(writer, false);
            }
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
