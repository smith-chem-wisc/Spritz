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
        public string ConfigFile { get; }
        public string ResourcesDirectory { get; }
        public string PathToWorkflow { get; }
        public string SpritzContainerName { get; set; }
        public string SnakemakeCommand { get; private set; }
        public string SpritzCMDCommand { get; set; }

        public static readonly string CurrentVersion = "0.3.0"; // should be the same here, in config.yaml, and in common.smk
        public static readonly bool PrebuiltSpritzMods = true; // always using prebuilt library now
        public RunnerEngine(Tuple<string, Options> task, string outputFolder)
        {
            // set up directories to mount to docker container as volumes
            AnalysisDirectory = task != null ? task.Item2.AnalysisDirectory : outputFolder;

            var pathToConfig = Path.Combine(AnalysisDirectory, "config");
            if (!Directory.Exists(pathToConfig))
            {
                Directory.CreateDirectory(pathToConfig);
            }
            ConfigDirectory = pathToConfig;
            ConfigFile = Path.Combine(ConfigDirectory, "config.yaml");

            var resourcesPath = Path.Combine(Path.GetDirectoryName(AnalysisDirectory), "resources");
            if (!Directory.Exists(resourcesPath))
            {
                Directory.CreateDirectory(resourcesPath);
            }
            ResourcesDirectory = resourcesPath;

            // path to workflow.txt
            PathToWorkflow = Path.Combine(AnalysisDirectory, "workflow_" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss") + ".txt");
            SpritzContainerName = $"spritz{PathToWorkflow.GetHashCode()}";
        }

        public string GenerateCommandsDry(string dockerImageName, string spritzCmdCommand)
        {
            string command = "";
            if (dockerImageName.Contains("smithlab"))
            { 
                command += $"docker pull {dockerImageName} ;";
            }
            command += $"docker run --rm -i -t --user=root --name {SpritzContainerName} " +
                $"-v \"\"\"{AnalysisDirectory}:/app/analysis" + "\"\"\" " +
                $"-v \"\"\"{ResourcesDirectory}:/app/resources" + "\"\"\" " +
                $"{dockerImageName} {spritzCmdCommand}; docker stop spritz{PathToWorkflow.GetHashCode()}";
            return command;
        }

        public string GenerateSpritzCMDCommand(Options options)
        {
            string command = $"/opt/conda/lib/dotnet/dotnet SpritzCMD.dll {SpritzCmdAppArgInfoStrings.GenerateSpritzCMDArgs(options)}";
            SpritzCMDCommand = command;
            return command;
        }

        public string GenerateSnakemakeCommand(Options options, bool setup)
        {
            string cmd = "";
            cmd += $"snakemake -j {options.Threads} --use-conda --conda-frontend mamba --configfile {Path.Combine(ConfigDirectory, "config.yaml")}";
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

        public void WriteConfig(Options options, string analysisDirectoryStr)
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
            if (options.Quantify) analysisStrings.Add("quant");

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
            analysisDirectory.Add(analysisDirectoryStr);
            rootMappingNode.Add("analysisDirectory", analysisDirectory);

            // process reference string
            var reference = options.Reference.Split(',');
            if (reference.Length != 4)
            {
                throw new SpritzException($"Error: the reference string \"{reference}\" does not have four comma-separated elements corresponding to a line from genomes.csv.");
            }
            string releaseStr = reference[0];
            string speciesStr = reference[1];
            string organismStr = reference[2];
            string referenceStr = reference[3];

            // write ensembl release
            YamlScalarNode release = new(releaseStr[8..]);
            release.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("release", release);

            // write species
            YamlScalarNode species = new(speciesStr.First().ToString().ToUpper() + speciesStr[1..]);
            species.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("species", species);

            // write organism
            YamlScalarNode organism = new(organismStr);
            organism.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("organism", organism);

            // write genome [e.g. GRCm38]
            YamlScalarNode genome = new(referenceStr);
            genome.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("genome", genome);

            // list the analyses to perform
            var analyses = new YamlSequenceNode();
            rootMappingNode.Add("analyses", AddParam(analysisStrings.ToArray(), analyses));

            // record the version of spritz
            YamlScalarNode version = new(RunnerEngine.CurrentVersion);
            version.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("spritzversion", version);

            // record that the spritzmods dll will be prebuilt via SpritzCMD
            rootMappingNode.Add("prebuilt_spritz_mods", new YamlScalarNode(PrebuiltSpritzMods ? "True" : "False"));

            using TextWriter writer = File.CreateText(ConfigFile);
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