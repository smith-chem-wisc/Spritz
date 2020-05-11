using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using YamlDotNet.Core;
using YamlDotNet.Core.Events;
using YamlDotNet.RepresentationModel;

namespace SpritzGUI
{
    public class EverythingRunnerEngine
    {
        private readonly Tuple<string, Options> task;
        private string outputFolder;

        public EverythingRunnerEngine(Tuple<string, Options> task, string outputFolder)
        {
            this.task = task;
            this.outputFolder = outputFolder;
        }

        public string Arguments { get; set; }
        public string StdErr { get; set; }
        public static string AnalysisDirectory { get; set; }
        public static string ConfigDirectory { get; set; }
        public static string DataDirectory { get; set; }
        public static string PathToWorkflow { get; set; }

        public void Run()
        {
            WriteConfig(task.Item2);

            Process proc = new Process();
            proc.StartInfo.FileName = "Powershell.exe";
            proc.StartInfo.Arguments = "docker pull smithlab/spritz ; docker run --rm -t -i --name spritz " +
                "-v \"\"\"" + AnalysisDirectory + ":/app/analysis" + "\"\"\" " +
                "-v \"\"\"" + DataDirectory + ":/app/data" + "\"\"\" " +
                "-v \"\"\"" + ConfigDirectory + ":/app/configs\"\"\" " +
                "smithlab/spritz > " + "\"\"\"" + PathToWorkflow + "\"\"\"";

            //proc.StartInfo.CreateNoWindow = true;
            proc.StartInfo.UseShellExecute = true;
            //proc.StartInfo.RedirectStandardError = true;
            proc.Start();
            proc.WaitForExit();
            //StdErr = proc.StandardError.ReadToEnd();
        }

        public IEnumerable<string> GenerateCommandsDry()
        {
            yield return "docker pull smithlab/spritz ; docker run --rm -t -i --name spritz " +
                "-v \"\"\"" + AnalysisDirectory + ":/app/analysis" + "\"\"\" " +
                "-v \"\"\"" + DataDirectory + ":/app/data" + "\"\"\" " +
                "-v \"\"\"" + ConfigDirectory + ":/app/configs\"\"\" " +
                "smithlab/spritz > " + "\"\"\"" + PathToWorkflow + "\"\"\"";
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

        public void SetUpDirectories()
        {
            // set up directories to mount to docker container as volumes
            AnalysisDirectory = task.Item2.AnalysisDirectory;

            var pathToConfig = Path.Combine(Directory.GetCurrentDirectory(), "configs");
            if (!Directory.Exists(pathToConfig))
            {
                Directory.CreateDirectory(pathToConfig);
            }
            ConfigDirectory = pathToConfig;

            var pathToDataFiles = Path.Combine(AnalysisDirectory, "data");
            if (!Directory.Exists(pathToDataFiles))
            {
                Directory.CreateDirectory(pathToDataFiles);
            }
            DataDirectory = pathToDataFiles;

            // path to workflow.txt
            PathToWorkflow = Path.Combine(AnalysisDirectory, "workflow.txt");
        }

        private void WriteConfig(Options options)
        {
            const string initialContent = "---\nversion: 1\n"; // needed to start writing yaml file

            var sr = new StringReader(initialContent);
            var stream = new YamlStream();
            stream.Load(sr);

            var rootMappingNode = (YamlMappingNode)stream.Documents[0].RootNode;

            var sras = options.SraAccession.Split(',');
            var fqs = options.Fastq1.Split(',') ?? new string[0];

            // write user input sras
            var accession = new YamlSequenceNode();
            rootMappingNode.Add("sra", AddParam(sras, accession));

            // write user defined analysis directory (input and output folder)
            var analysisDirectory = new YamlSequenceNode();
            analysisDirectory.Style = SequenceStyle.Flow;
            analysisDirectory.Add("analysis");
            rootMappingNode.Add("analysisDirectory", analysisDirectory);

            // write user input fastqs
            var fq = new YamlSequenceNode();
            rootMappingNode.Add("fq", AddParam(fqs, fq));

            // write ensembl release
            var release = new YamlScalarNode(options.Release.Substring(8));
            release.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("release", release);

            // write species
            var species = new YamlScalarNode(options.Species.First().ToString().ToUpper() + options.Species.Substring(1));
            species.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("species", species);

            // write species
            var organism = new YamlScalarNode(options.Organism);
            organism.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("organism", organism);

            // write genome [e.g. GRCm38]
            var genome = new YamlScalarNode(options.Reference);
            genome.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("genome", genome);

            // write snpeff (hardcoded for now)
            var snpeff = new YamlScalarNode(options.SnpEff);
            snpeff.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("snpeff", snpeff);

            using (TextWriter writer = File.CreateText(Path.Combine(ConfigDirectory, "config.yaml")))
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