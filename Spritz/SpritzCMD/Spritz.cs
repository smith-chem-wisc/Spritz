﻿using Fclp;
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
        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to Spritz!");
            FluentCommandLineParser<SpritzCmdAppArguments> p = new();

            Options defaults = new(Environment.ProcessorCount);
            p.Setup(arg => arg.AnalysisDirectory)
                .As(SpritzCmdAppArgInfoStrings.AnalysisDirectoryShort,
                    SpritzCmdAppArgInfoStrings.AnalysisDirectoryLong)
                .SetDefault(defaults.AnalysisDirectory)
                .WithDescription(SpritzCmdAppArgInfoStrings.AnalysisDirectoryDesc);

            p.Setup(arg => arg.AnalyzeVariants)
                .As(SpritzCmdAppArgInfoStrings.AnalyzeVariantsShort,
                    SpritzCmdAppArgInfoStrings.AnalyzeVariantsLong)
                .SetDefault(defaults.AnalyzeVariants)
                .WithDescription(SpritzCmdAppArgInfoStrings.AnalyzeVariantsDesc);

            p.Setup(arg => arg.AnalyzeIsoforms)
                .As(SpritzCmdAppArgInfoStrings.AnalyzeIsoformsShort,
                    SpritzCmdAppArgInfoStrings.AnalyzeIsoformsLong)
                .SetDefault(defaults.AnalyzeIsoforms)
                .WithDescription(SpritzCmdAppArgInfoStrings.AnalyzeIsoformsDesc);

            p.Setup(arg => arg.Quantify)
                .As(SpritzCmdAppArgInfoStrings.QuantifyShort,
                    SpritzCmdAppArgInfoStrings.QuantifyLong)
                .SetDefault(defaults.Quantify)
                .WithDescription(SpritzCmdAppArgInfoStrings.QuantifyDesc);

            p.Setup(arg => arg.AvailableReferences)
                .As(SpritzCmdAppArgInfoStrings.AvailableReferencesShort,
                    SpritzCmdAppArgInfoStrings.AvailableReferencesLong)
                .SetDefault(false)
                .WithDescription(SpritzCmdAppArgInfoStrings.AvailableReferencesDesc);

            p.Setup(arg => arg.AnalysisSetup)
                .As(SpritzCmdAppArgInfoStrings.AnalysisSetupShort,
                    SpritzCmdAppArgInfoStrings.AnalysisSetupLong)
                .SetDefault(false)
                .WithDescription(SpritzCmdAppArgInfoStrings.AnalysisSetupDesc);

            p.Setup(arg => arg.Fastq1)
                .As(SpritzCmdAppArgInfoStrings.Fastq1Short,
                    SpritzCmdAppArgInfoStrings.Fastq1Long)
                .WithDescription(SpritzCmdAppArgInfoStrings.Fastq1Desc);

            p.Setup(arg => arg.Fastq2)
                .As(SpritzCmdAppArgInfoStrings.Fastq2Short,
                    SpritzCmdAppArgInfoStrings.Fastq2Long)
                .WithDescription(SpritzCmdAppArgInfoStrings.Fastq2Desc);

            p.Setup(arg => arg.Fastq1SingleEnd)
                .As(SpritzCmdAppArgInfoStrings.Fastq1SingleEndShort,
                    SpritzCmdAppArgInfoStrings.Fastq1SingleEndLong)
                .WithDescription(SpritzCmdAppArgInfoStrings.Fastq1SingleEndDesc);

            p.Setup(arg => arg.SraAccession)
                .As(SpritzCmdAppArgInfoStrings.SraAccessionShort,
                    SpritzCmdAppArgInfoStrings.SraAccessionLong)
                .WithDescription(SpritzCmdAppArgInfoStrings.SraAccessionDesc);

            p.Setup(arg => arg.SraAccessionSingleEnd)
                .As(SpritzCmdAppArgInfoStrings.SraAccessionSingleEndShort,
                    SpritzCmdAppArgInfoStrings.SraAccessionSingleEndLong)
                .WithDescription(SpritzCmdAppArgInfoStrings.SraAccessionSingleEndDesc);

            p.Setup(arg => arg.Threads)
                .As(SpritzCmdAppArgInfoStrings.ThreadsShort,
                    SpritzCmdAppArgInfoStrings.ThreadsLong)
                .SetDefault(defaults.Threads)
                .WithDescription(SpritzCmdAppArgInfoStrings.ThreadsDesc);

            p.Setup(arg => arg.Reference)
                 .As(SpritzCmdAppArgInfoStrings.ReferenceShort,
                    SpritzCmdAppArgInfoStrings.ReferenceLong)
                .WithDescription(SpritzCmdAppArgInfoStrings.ReferenceDesc);

            string helpoutro = "";
            helpoutro += $"The Spritz commandline interface intended to be run within a conda environment containing the programs snakemake and mamba." + Environment.NewLine;
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
            var possibleMatches = typeof(SpritzCmdAppArguments).GetFields(BindingFlags.Instance | BindingFlags.Static | BindingFlags.NonPublic | BindingFlags.NonPublic);
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
            else if (p.Object.AvailableReferences)
            {
                Console.WriteLine();
                Console.WriteLine($"Saving the list of available references to {Path.Combine(analysisDirectory, "genomes.csv")}.");
                string genomesPath = Path.Combine(Directory.GetCurrentDirectory(), "genomes.csv");
                Directory.CreateDirectory(analysisDirectory);
                string dest = Path.Combine(analysisDirectory, Path.GetFileName(genomesPath));
                if (File.Exists(dest))
                {
                    Console.WriteLine($"File {dest} already exists. Please check it out there.");
                }
                else
                {
                    File.Copy(genomesPath, dest);
                }
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

                Options options = ParseOptions(p.Object, analysisDirectory);

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
            }
        }

        private static Options ParseOptions(SpritzCmdAppArguments aa, string analysisDirectory)
        {
            Options options = new(aa.Threads);
            options.AnalysisDirectory = analysisDirectory;
            options.Fastq1 = aa.Fastq1 ?? "";
            options.Fastq2 = aa.Fastq2 ?? "";
            options.Fastq1SingleEnd = aa.Fastq1SingleEnd ?? "";
            options.SraAccession = aa.SraAccession ?? "";
            options.SraAccessionSingleEnd = aa.SraAccessionSingleEnd ?? "";
            options.Threads = aa.Threads;
            options.Reference = aa.Reference;
            options.AnalyzeVariants = aa.AnalyzeVariants;
            options.AnalyzeIsoforms = aa.AnalyzeIsoforms;
            options.Quantify = aa.Quantify;
            return options;
        }
    }
}