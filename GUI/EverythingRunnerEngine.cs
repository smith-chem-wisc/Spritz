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
