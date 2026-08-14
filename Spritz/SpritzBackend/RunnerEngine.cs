using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
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

        /// <summary>
        /// The version this build reports, and the Docker Hub tag it pulls. Read from the assembly rather
        /// than hardcoded, so the MSBuild Version property is the only place it is set: Directory.Build.props
        /// holds the development default and release.yml overrides it with /p:Version from the git tag.
        /// Previously this literal had to be edited in step with Product.wxs and the tag by hand.
        /// </summary>
        public static readonly string CurrentVersion = ReadCompiledVersion();

        /// <summary>
        /// The Ensembl release number from a reference string's first element, which genomes.csv writes as
        /// "release-116". Accepts a bare number too.
        /// </summary>
        private static string EnsemblReleaseNumber(string releaseField, string wholeReference)
        {
            const string prefix = "release-";
            string number = releaseField.StartsWith(prefix, StringComparison.OrdinalIgnoreCase)
                ? releaseField[prefix.Length..]
                : releaseField;

            if (number.Length == 0 || !number.All(char.IsDigit))
            {
                throw new SpritzException($"Error: \"{releaseField}\" in the reference string " +
                    $"\"{wholeReference}\" is not an Ensembl release. Expected something like " +
                    $"\"release-116\" or \"116\".");
            }

            return number;
        }

        private static string ReadCompiledVersion()
        {
            // InformationalVersion carries the full Version string. Source-link style builds append "+<sha>",
            // which is not part of the version and is not in the Docker tag, so it is trimmed.
            string informational = typeof(RunnerEngine).Assembly
                .GetCustomAttribute<AssemblyInformationalVersionAttribute>()?.InformationalVersion;
            if (!string.IsNullOrWhiteSpace(informational))
            {
                return informational.Split('+')[0];
            }

            // AssemblyVersion is always present, but is four-part, so take the first three to match the tag.
            return typeof(RunnerEngine).Assembly.GetName().Version?.ToString(3) ?? "0.0.0";
        }

        public static readonly bool PrebuiltSpritzMods = true; // always using prebuilt library now
        public RunnerEngine(Tuple<string, SpritzOptions> task, string outputFolder)
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
            string runStamp = DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss");
            PathToWorkflow = Path.Combine(AnalysisDirectory, "workflow_" + runStamp + ".txt");

            // Derived from the same timestamp as the workflow file, so the two correspond and the name
            // is readable in `docker ps`. String.GetHashCode was used here, which is not stable across
            // processes in .NET Core - a second process could not name the container this one started.
            SpritzContainerName = $"spritz-{runStamp}";
        }

        /// <summary>Registries Spritz publishes to. An image from one of these is pulled before running.</summary>
        private static readonly string[] PublishedImagePrefixes =
        {
            "smithlab/", "docker.io/smithlab/", "ghcr.io/smith-chem-wisc/",
        };

        /// <summary>Set to pull an image that is not on a known published registry, e.g. on a fork.</summary>
        public bool AlwaysPull { get; set; }

        /// <summary>
        /// Whether to pull before running. Matched on an exact prefix rather than a substring: a local
        /// image called "smithlab-test" used to match and be overwritten by the published one, and a
        /// fork publishing anywhere else was never pulled at all.
        /// </summary>
        public bool ShouldPull(string imageName) =>
            AlwaysPull || PublishedImagePrefixes.Any(prefix =>
                imageName.StartsWith(prefix, StringComparison.OrdinalIgnoreCase));

        /// <summary>The container runtime to drive. Podman by default; see ContainerRuntime.</summary>
        public ContainerRuntime Runtime { get; set; } = ContainerRuntime.Podman;

        public string GenerateCommandsDry(string dockerImageName, string spritzCmdCommand)
        {
            string imageWithVersion = dockerImageName.Contains(":") ? dockerImageName : $"{dockerImageName}:{CurrentVersion}";
            return Runtime == ContainerRuntime.Apptainer
                ? ApptainerCommand(imageWithVersion, spritzCmdCommand)
                : DockerCompatibleCommand(imageWithVersion, spritzCmdCommand);
        }

        /// <summary>
        /// The run command as an executable plus a list of arguments, which is what
        /// ProcessStartInfo.ArgumentList takes. Nothing is quoted or escaped here, so a path containing
        /// a quote or a space cannot change what the command means.
        ///
        /// GenerateCommandsDry below produces the same thing as a single string for the WPF GUI, which
        /// shells through Powershell.exe and therefore has to escape. Prefer this in new callers; the
        /// string form can go when that GUI does.
        /// </summary>
        public (string Executable, IReadOnlyList<string> Arguments) GenerateRunArguments(string imageName, string spritzCmdCommand)
        {
            string imageWithVersion = imageName.Contains(":") ? imageName : $"{imageName}:{CurrentVersion}";
            var arguments = new List<string>();

            if (Runtime == ContainerRuntime.Apptainer)
            {
                arguments.AddRange(new[] { "run", "--cleanenv", "--writable-tmpfs", "--pwd", "/app/spritz/" });
                arguments.AddRange(new[] { "--bind", $"{AnalysisDirectory}:/app/spritz/results/" });
                arguments.AddRange(new[] { "--bind", $"{ResourcesDirectory}:/app/spritz/resources" });
                arguments.Add(imageWithVersion.StartsWith("docker://") ? imageWithVersion : $"docker://{imageWithVersion}");
            }
            else
            {
                arguments.AddRange(new[] { "run", "--rm", "-i", "-t", "--user=root", "--name", SpritzContainerName });
                arguments.AddRange(new[] { "-v", $"{AnalysisDirectory}:/app/spritz/results/" });
                arguments.AddRange(new[] { "-v", $"{ResourcesDirectory}:/app/spritz/resources" });
                arguments.Add(imageWithVersion);
            }

            arguments.AddRange(spritzCmdCommand.Split(' ', StringSplitOptions.RemoveEmptyEntries));
            return (Runtime.Executable(), arguments);
        }

        /// <summary>
        /// Docker and Podman take the same flags for everything used here, so one builder serves both.
        /// </summary>
        private string DockerCompatibleCommand(string imageWithVersion, string spritzCmdCommand)
        {
            string exe = Runtime.Executable();
            string command = ShouldPull(imageWithVersion) ? $"{exe} pull {imageWithVersion};" : "";
            command +=
                $"{exe} run --rm -i -t --user=root --name {SpritzContainerName} " +
                $"-v \"\"\"{AnalysisDirectory}:/app/spritz/results/\"\"\" " +
                $"-v \"\"\"{ResourcesDirectory}:/app/spritz/resources\"\"\" " +
                $"{imageWithVersion} {spritzCmdCommand}; " +
                $"{exe} stop {SpritzContainerName}";
            return command;
        }

        /// <summary>
        /// Apptainer has no daemon, so there is no container to name or stop, and it binds rather than
        /// mounts. It reads an OCI image through a docker:// URI, so the same published image is used.
        ///
        /// --writable-tmpfs is required because the image is read-only and the workflow writes inside
        /// it. That tmpfs is memory-backed, so a long run wants a real directory bound over the
        /// snakemake state instead - see the cluster guide.
        /// </summary>
        private string ApptainerCommand(string imageWithVersion, string spritzCmdCommand)
        {
            string uri = imageWithVersion.StartsWith("docker://") ? imageWithVersion : $"docker://{imageWithVersion}";
            return
                // --pwd because Apptainer ignores the image WORKDIR and starts in the host directory,
                // so SpritzCMD.dll would not be found.
                $"apptainer run --cleanenv --writable-tmpfs --pwd /app/spritz/ " +
                $"--bind \"\"\"{AnalysisDirectory}:/app/spritz/results/\"\"\" " +
                $"--bind \"\"\"{ResourcesDirectory}:/app/spritz/resources\"\"\" " +
                $"{uri} {spritzCmdCommand}";
        }

        public string GenerateSpritzCMDCommand(SpritzOptions options)
        {
            string command = $"conda run --no-capture-output --live-stream dotnet SpritzCMD.dll {SpritzOptionStrings.GenerateSpritzCMDArgs(options)}";
            SpritzCMDCommand = command;
            return command;
        }

        public string GenerateSnakemakeCommand(SpritzOptions options, bool setup)
        {
            string cmd = "";
            // No --conda-frontend: snakemake 9 accepts the flag but prints "Ignoring the alternative
            // conda frontend setting (mamba)" and uses conda, which now solves via libmamba anyway.
            // Passing it only produced a warning on every run.
            cmd += $"snakemake -j {options.Threads} --use-conda --configfile {Path.Combine(ConfigDirectory, "config.yaml")}";
            if (setup)
            {
                // The rule is named "setup" and writes ../resources/setup.txt. "setup.txt" matches
                // neither, so snakemake refused with MissingRuleException.
                cmd += " setup";
            }
            SnakemakeCommand = cmd;
            return cmd;
        }

        public string GenerateTopComand()
        {
            // Apptainer has no daemon, so there is nothing to inspect this way.
            return Runtime.IsDockerCompatible()
                ? $"{Runtime.Executable()} container top {SpritzContainerName}"
                : string.Empty;
        }

        public void WriteConfig(SpritzOptions options, string analysisDirectoryStr)
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
                throw new SpritzException($"Error: the reference string \"{options.Reference}\" has " +
                    $"{reference.Length} comma-separated element(s), not four. Copy a whole line from " +
                    $"genomes.csv, e.g. \"release-116,homo_sapiens,human,GRCh38\".");
            }
            string releaseStr = reference[0];
            string speciesStr = reference[1];
            string organismStr = reference[2];
            string referenceStr = reference[3];

            for (int i = 0; i < reference.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(reference[i]))
                {
                    throw new SpritzException($"Error: element {i + 1} of the reference string " +
                        $"\"{options.Reference}\" is empty. Copy a whole line from genomes.csv.");
                }
            }

            // write ensembl release
            YamlScalarNode release = new(EnsemblReleaseNumber(releaseStr, options.Reference));
            release.Style = ScalarStyle.DoubleQuoted;
            rootMappingNode.Add("release", release);

            // write species; invariant casing, and the field is guarded non-empty above
            YamlScalarNode species = new(char.ToUpperInvariant(speciesStr[0]) + speciesStr[1..]);
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
                string testDirectory = Path.Combine(path, $"TestSpritzPermissions{Guid.NewGuid():N}");
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
