using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using YamlDotNet.Core;
using YamlDotNet.Core.Events;
using YamlDotNet.RepresentationModel;

namespace SpritzBackend
{
    public class RunnerEngine
    {
        public string Arguments { get; set; }
        public string AnalysisDirectory { get; }
        public string ConfigDirectory { get; }
        public string DataDirectory { get; }
        public string PathToWorkflow { get; }
        public string SpritzContainerName { get; set; }
        public string SnakemakeCommand { get; private set; }

        public RunnerEngine(Tuple<string, Options> task, string outputFolder)
        {
            // set up directories to mount to docker container as volumes
            AnalysisDirectory = task != null ? task.Item2.AnalysisDirectory : outputFolder;

            var pathToConfig = Path.Combine(AnalysisDirectory, "configs");
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
            PathToWorkflow = Path.Combine(AnalysisDirectory, "workflow_" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss") + ".txt");
            SpritzContainerName = $"spritz{PathToWorkflow.GetHashCode()}";
        }

        public string GenerateCommandsDry(string dockerImageName, string snakemakeCommand)
        {
            string command = "";
            if (dockerImageName.Contains("smithlab"))
            { 
                command += $"docker pull {dockerImageName} ;";
            }
            command += $"docker run --rm -i -t --user=root --name {SpritzContainerName} " +
                $"-v \"\"\"{AnalysisDirectory}:/app/analysis" + "\"\"\" " +
                $"-v \"\"\"{DataDirectory}:/app/resources" + "\"\"\" " +
                $"-v \"\"\"{ConfigDirectory}:/app/configs\"\"\" " +
                $"{dockerImageName} {snakemakeCommand}; docker stop spritz{PathToWorkflow.GetHashCode()}";
            return command;
        }

        public string GenerateSnakemakeCommand(Options options, bool setup)
        {
            string cmd = "";
            cmd += $"snakemake -j {options.Threads} --use-conda --conda-frontend mamba";
            if (setup)
            {
                cmd += " setup.txt";
            }
            SnakemakeCommand = cmd;
            return cmd;
        }

        public string GenerateTopComand()
        {
            return $"docker container top spritz{PathToWorkflow.GetHashCode()}";
        }

        public void WriteConfig(Options options)
        {
            const string initialContent = "---\nversion: 1\n"; // needed to start writing yaml file

            var sr = new StringReader(initialContent);
            var stream = new YamlStream();
            stream.Load(sr);

            var rootMappingNode = (YamlMappingNode)stream.Documents[0].RootNode;

            var sras = options.SraAccession.Split(',');
            var sras_se = options.SraAccessionSingleEnd.Split(',');
            var fqs = options.Fastq1.Split(',') ?? Array.Empty<string>();
            var fqs_se = options.Fastq1SingleEnd.Split(',') ?? Array.Empty<string>();
            var analysisStrings = new List<string>();
            if (options.AnalyzeVariants) analysisStrings.Add("variant");
            if (options.AnalyzeIsoforms) analysisStrings.Add("isoform");

            // write user input paired-end sras
            YamlSequenceNode accession = new();
            rootMappingNode.Add("sra", AddParam(sras, accession));

            // write user input paired-end sras
            YamlSequenceNode accession_se = new();
            rootMappingNode.Add("sra_se", AddParam(sras_se, accession_se));

            // write user input paired-end fastqs
            YamlSequenceNode fq = new();
            rootMappingNode.Add("fq", AddParam(fqs, fq));

            // write user input paired-end fastqs
            YamlSequenceNode fq_se = new();
            rootMappingNode.Add("fq_se", AddParam(fqs_se, fq_se));

            // write user defined analysis directory (input and output folder)
            YamlSequenceNode analysisDirectory = new();
            analysisDirectory.Style = SequenceStyle.Flow;
            analysisDirectory.Add("analysis");
            rootMappingNode.Add("analysisDirectory", analysisDirectory);

            // write ensembl release
            YamlScalarNode release = new(options.Release[8..]);
            release.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("release", release);

            // write species
            YamlScalarNode species = new(options.Species.First().ToString().ToUpper() + options.Species[1..]);
            species.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("species", species);

            // write species
            YamlScalarNode organism = new(options.Organism);
            organism.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("organism", organism);

            // write genome [e.g. GRCm38]
            YamlScalarNode genome = new(options.Reference);
            genome.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("genome", genome);

            // list the analyses to perform
            var analyses = new YamlSequenceNode();
            rootMappingNode.Add("analyses", AddParam(analysisStrings.ToArray(), analyses));

            // record the version of spritz
            YamlScalarNode version = new(options.SpritzVersion);
            version.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("spritzversion", version);

            using TextWriter writer = File.CreateText(Path.Combine(ConfigDirectory, "config.yaml"));
            stream.Save(writer, false);
        }

        private static YamlSequenceNode AddParam(string[] items, YamlSequenceNode node)
        {
            node.Style = SequenceStyle.Flow;
            foreach (string item in items.Where(x => x.Length > 0))
            {
                node.Add(item);
            }
            return node;
        }

        public static bool IsDirectoryWritable(string path)
        {
            try
            {
                string testDirectory = Path.Combine(path, $"TestSpritzPermissions{path.GetHashCode()}");
                Directory.CreateDirectory(testDirectory);
                Directory.Delete(testDirectory);
            }
            catch (Exception)
            {
                return false;
            }
            return true;
        }

        public static string TrimQuotesOrNull(string a)
        {
            return a == null ? a : a.Trim('"');
        }
    }
}