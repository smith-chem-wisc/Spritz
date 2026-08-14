using Fclp;
using SpritzBackend;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;

namespace SpritzCMD
{
    internal class Spritz
    {
        /// <summary>
        /// Data lines in genomes.csv, skipping the generator's comment header.
        /// </summary>
        private static int CountReferences(string path)
        {
            int count = 0;
            foreach (string line in File.ReadLines(path))
            {
                if (line.Length > 0 && !line.StartsWith("#", StringComparison.Ordinal))
                {
                    count++;
                }
            }
            return count;
        }

        /// <summary>
        /// Stops the parser reading "/" as a Windows-style option prefix, which makes absolute POSIX paths
        /// unusable as option values. Replaced with a duplicate of an existing prefix, not a sentinel: a
        /// one-character prefix matches every argument and strips its first character.
        /// </summary>
        /// <summary>
        /// Rewrites the installed genomes.csv from the Ensembl FTP listings, by running the generator that
        /// produces the checked-in copy. Existing rows are kept, so reference strings already in use keep
        /// working. Needs network access and python, which the conda environment Spritz documents provides.
        /// </summary>
        private static void RefreshGenomesFromEnsembl(string destination)
        {
            string script = Path.Combine(AppContext.BaseDirectory, "workflow", "scripts", "update_genomes.py");
            if (!File.Exists(script))
            {
                throw new SpritzException(
                    $"Error: cannot refresh the reference list because '{script}' is missing from the " +
                    $"Spritz installation. Use {SpritzOptionStrings.AvailableReferencesShort} for the list " +
                    $"that shipped with this version.");
            }

            Console.WriteLine($"Refreshing {destination} from Ensembl. This reads the FTP listings for " +
                $"every species and takes a few minutes; progress follows.");

            // python3 first: a Linux box commonly has python3 and no python at all, while a conda
            // environment has both. Output is not redirected, so the generator's progress reaches the
            // console as it happens rather than arriving after several silent minutes.
            foreach (string interpreter in new[] { "python3", "python" })
            {
                var start = new ProcessStartInfo(interpreter) { UseShellExecute = false };
                start.ArgumentList.Add(script);
                start.ArgumentList.Add("--output");
                start.ArgumentList.Add(destination);

                Process process;
                try
                {
                    process = Process.Start(start);
                }
                catch (System.ComponentModel.Win32Exception)
                {
                    continue; // not on the PATH under this name
                }

                using (process)
                {
                    process.WaitForExit();
                    if (process.ExitCode != 0)
                    {
                        throw new SpritzException(
                            $"Error: refreshing the reference list from Ensembl failed (exit code " +
                            $"{process.ExitCode}); see the output above. The list that shipped with this " +
                            $"version is still available with -{SpritzOptionStrings.AvailableReferencesShort}.");
                    }
                }

                return;
            }

            throw new SpritzException(
                $"Error: neither python3 nor python is on the PATH, so the reference list cannot be " +
                $"refreshed. Run Spritz from an environment that has python, or use " +
                $"-{SpritzOptionStrings.AvailableReferencesShort} for the list that shipped with it.");
        }

        private static void AllowAbsolutePosixPaths()
        {
            string[] prefixes = Fclp.Internals.SpecialCharacters.OptionPrefix;
            for (int i = 0; i < prefixes.Length; i++)
            {
                if (prefixes[i] == "/")
                {
                    prefixes[i] = "--";
                }
            }
        }

        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to Spritz!");
            AllowAbsolutePosixPaths();
            FluentCommandLineParser<SpritzOptions> p = new();

            // Get defaults
            SpritzOptions defaults = new();
            defaults.AnalysisDirectory = SpritzOptions.DefaultAnalysisDirectory();
            defaults.Threads = Environment.ProcessorCount;

            p.Setup(arg => arg.AnalysisDirectory)
                .As(SpritzOptionStrings.AnalysisDirectoryShort,
                    SpritzOptionStrings.AnalysisDirectoryLong)
                .SetDefault(defaults.AnalysisDirectory)
                .WithDescription(SpritzOptionStrings.AnalysisDirectoryDesc);

            p.Setup(arg => arg.AnalyzeVariants)
                .As(SpritzOptionStrings.AnalyzeVariantsShort,
                    SpritzOptionStrings.AnalyzeVariantsLong)
                .SetDefault(defaults.AnalyzeVariants)
                .WithDescription(SpritzOptionStrings.AnalyzeVariantsDesc);

            p.Setup(arg => arg.AnalyzeIsoforms)
                .As(SpritzOptionStrings.AnalyzeIsoformsShort,
                    SpritzOptionStrings.AnalyzeIsoformsLong)
                .SetDefault(defaults.AnalyzeIsoforms)
                .WithDescription(SpritzOptionStrings.AnalyzeIsoformsDesc);

            p.Setup(arg => arg.Quantify)
                .As(SpritzOptionStrings.QuantifyShort,
                    SpritzOptionStrings.QuantifyLong)
                .SetDefault(defaults.Quantify)
                .WithDescription(SpritzOptionStrings.QuantifyDesc);

            p.Setup(arg => arg.AvailableReferences)
                .As(SpritzOptionStrings.AvailableReferencesShort,
                    SpritzOptionStrings.AvailableReferencesLong)
                .SetDefault(false)
                .WithDescription(SpritzOptionStrings.AvailableReferencesDesc);

            p.Setup(arg => arg.FetchGenomes)
                .As(SpritzOptionStrings.FetchGenomesShort,
                    SpritzOptionStrings.FetchGenomesLong)
                .SetDefault(false)
                .WithDescription(SpritzOptionStrings.FetchGenomesDesc);

            p.Setup(arg => arg.AnalysisSetup)
                .As(SpritzOptionStrings.AnalysisSetupShort,
                    SpritzOptionStrings.AnalysisSetupLong)
                .SetDefault(false)
                .WithDescription(SpritzOptionStrings.AnalysisSetupDesc);

            p.Setup(arg => arg.Fastq1)
                .As(SpritzOptionStrings.Fastq1Short,
                    SpritzOptionStrings.Fastq1Long)
                .WithDescription(SpritzOptionStrings.Fastq1Desc);

            p.Setup(arg => arg.Fastq2)
                .As(SpritzOptionStrings.Fastq2Short,
                    SpritzOptionStrings.Fastq2Long)
                .WithDescription(SpritzOptionStrings.Fastq2Desc);

            p.Setup(arg => arg.Fastq1SingleEnd)
                .As(SpritzOptionStrings.Fastq1SingleEndShort,
                    SpritzOptionStrings.Fastq1SingleEndLong)
                .WithDescription(SpritzOptionStrings.Fastq1SingleEndDesc);

            p.Setup(arg => arg.SraAccession)
                .As(SpritzOptionStrings.SraAccessionShort,
                    SpritzOptionStrings.SraAccessionLong)
                .WithDescription(SpritzOptionStrings.SraAccessionDesc);

            p.Setup(arg => arg.SraAccessionSingleEnd)
                .As(SpritzOptionStrings.SraAccessionSingleEndShort,
                    SpritzOptionStrings.SraAccessionSingleEndLong)
                .WithDescription(SpritzOptionStrings.SraAccessionSingleEndDesc);

            p.Setup(arg => arg.Threads)
                .As(SpritzOptionStrings.ThreadsShort,
                    SpritzOptionStrings.ThreadsLong)
                .SetDefault(defaults.Threads)
                .WithDescription(SpritzOptionStrings.ThreadsDesc);

            p.Setup(arg => arg.ContainerRuntime)
                .As(SpritzOptionStrings.ContainerRuntimeLong)
                .SetDefault("podman")
                .WithDescription(SpritzOptionStrings.ContainerRuntimeDesc);

            p.Setup(arg => arg.Reference)
                 .As(SpritzOptionStrings.ReferenceShort,
                    SpritzOptionStrings.ReferenceLong)
                .WithDescription(SpritzOptionStrings.ReferenceDesc);

            string helpoutro = "";
            helpoutro += $"The Spritz commandline interface intended to be run within a conda environment containing the programs snakemake and conda." + Environment.NewLine;
            helpoutro += Environment.NewLine;
            helpoutro += $"Example workflow using this tool:" + Environment.NewLine;
            helpoutro += $"1) Check out the available references with the -x command. Specify a target directory with -a." + Environment.NewLine;
            helpoutro += $"2) Run spritz with -r based on the genomes.csv file generated by 1), " + Environment.NewLine;
            helpoutro += $"and choose the workflow options -b to analyze variants, -c to analyze isoforms, or both, and results will be saved at directory specified by -a. " + Environment.NewLine;
            helpoutro += $"2b) Alternatively, specify false for both -v and -w to generate a reference proteogenomic database from the Ensembl references." + Environment.NewLine;
            helpoutro += Environment.NewLine;
            p.SetupHelp("h", "help")
                .Callback(text => Console.WriteLine(text + helpoutro));

            var result = p.Parse(args);

            // handle unrecognized and unmatched
            bool anyUnrecognized = result.AdditionalOptionsFound.Any();
            int countUnmatched = result.UnMatchedOptions.Count();
            var possibleMatches = typeof(SpritzOptions).GetFields(BindingFlags.Instance | BindingFlags.Static | BindingFlags.NonPublic | BindingFlags.NonPublic);
            if (anyUnrecognized)
            {
                throw new SpritzException($"Error: unrecognized commandline argument(s): {string.Join(',', result.AdditionalOptionsFound.Select(x => x.ToString()))}");
            }
            else if (countUnmatched == possibleMatches.Length)
            {
                result = p.Parse(new[] { "-h" });
            }

            string analysisDirectory = RunnerEngine.TrimQuotesOrNull(p.Object.AnalysisDirectory);
            Console.WriteLine($"Testing analysis directory {analysisDirectory}");
            if (!RunnerEngine.IsDirectoryWritable(analysisDirectory))
            {
                analysisDirectory = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "spritz", "results");
            }
            Console.WriteLine($"Using analysis directory {analysisDirectory}");
            Console.WriteLine();

            bool noSequencesSpecified =
                p.Object.SraAccession == null && p.Object.SraAccessionSingleEnd == null &&
                p.Object.Fastq1 == null && p.Object.Fastq1SingleEnd == null && p.Object.Fastq2 == null;
            bool analysisSpecified =
                p.Object.AnalyzeVariants || p.Object.AnalyzeIsoforms || p.Object.Quantify;

            if (result.HelpCalled)
            {
                return;
            }
            else if (p.Object.AvailableReferences || p.Object.FetchGenomes)
            {
                Console.WriteLine();

                // beside the assembly, where SpritzBackend.csproj copies it, not in the working directory
                string genomesPath = ReferenceString.InstalledGenomesPath();
                if (!File.Exists(genomesPath))
                {
                    throw new SpritzException(
                        $"Error: the list of available references is missing from the Spritz installation " +
                        $"(expected '{genomesPath}'). Reinstalling should restore it.");
                }

                Directory.CreateDirectory(analysisDirectory);
                string dest = Path.Combine(analysisDirectory, Path.GetFileName(genomesPath));

                // overwrite: an older copy from a previous version is what the user asked to refresh
                bool replaced = File.Exists(dest);
                File.Copy(genomesPath, dest, overwrite: true);

                // Written into the analysis directory, never over the installed copy: an installation can
                // sit somewhere the user cannot write. The generator merges with what is already there, so
                // copying first keeps the rows that shipped and any older releases with them.
                if (p.Object.FetchGenomes)
                {
                    RefreshGenomesFromEnsembl(dest);
                }

                Console.WriteLine($"{(replaced ? "Replaced" : "Saved")} the list of available references at {dest}.");
                Console.WriteLine($"It lists {CountReferences(dest)} references. Pass one line back with the " +
                    $"{SpritzOptionStrings.ReferenceShort} flag.");
                return;
            }
            else if (p.Object.Reference == null)
            {
                throw new SpritzException("Error: No reference specified. Please specify one with the -r flag that has four elements corresponding to a line from genomes.csv.");
            }
            else if (analysisSpecified && noSequencesSpecified)
            {
                throw new SpritzException("Error: An analysis was specified, but no sequences were specified to analyze. Please try again after specifying fastqs or sras.");
            }
            else
            {
                if (!analysisSpecified && noSequencesSpecified)
                    Console.WriteLine("NB: No sequences or analyses were specified, and so a reference database will be generated from Ensembl references only.");

                SpritzOptions options = CleanOptions(p.Object, analysisDirectory);

                RunnerEngine runner = new(new("", options), analysisDirectory);
                runner.WriteConfig(options, analysisDirectory);
                runner.GenerateSnakemakeCommand(options, p.Object.AnalysisSetup);
                string snakemakeArguments = runner.SnakemakeCommand["snakemake ".Length..];
                Console.WriteLine($"Running `{runner.SnakemakeCommand}`.");

                Process proc = new();
                proc.StartInfo.FileName = "snakemake";
                proc.StartInfo.Arguments = snakemakeArguments;
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.CreateNoWindow = true;
                proc.StartInfo.WorkingDirectory = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "workflow");
                proc.Start();
                proc.WaitForExit();
                Environment.ExitCode = proc.ExitCode;
            }
        }

        private static SpritzOptions CleanOptions(SpritzOptions aa, string analysisDirectory)
        {
            aa.AnalysisDirectory = analysisDirectory;
            aa.Reference ??= "";
            aa.Fastq1 ??= "";
            aa.Fastq2 ??= "";
            aa.Fastq1SingleEnd ??= "";
            aa.SraAccession ??= "";
            aa.SraAccessionSingleEnd ??= "";
            return aa;
        }
    }
}