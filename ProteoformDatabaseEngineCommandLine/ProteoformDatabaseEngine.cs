using System;
using System.Collections.Generic;
using System.Linq;
using RNASeqAnalysisWrappers;
using CommandLine;
using System.IO;
using System.Reflection;

namespace ProteoformDatabaseEngineCommandLine
{
    class ProteoformDatabaseEngine
    {
        public static void Main(string[] args)
        {

            #region Setup

            if (args.Contains("setup"))
            {
                WrapperUtility.Install(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
                return;
            }

            Options options = new Options();
            bool isValid = Parser.Default.ParseArgumentsStrict(args, options);

            #endregion Setup

            #region STAR Fusion Testing

            // TODO make a star fusion 2 protein runner instead of this mess...
            bool validStarFusionTest = options.AnalysisDirectory != null
                && options.BinDirectory != null
                && (options.Fastq1 != null || options.st)
            if (options.Command == "starFusionTest" && !starFusionRequirements.Any(x => x == null))
            {
                Directory.CreateDirectory(Path.Combine(options.AnalysisDirectory, "fusion_out"));
                STARFusionWrapper.Install(options.BinDirectory);
                STARFusionWrapper.RunStarFusion(options.BinDirectory,
                    "grch37",
                    8,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR791578_hg19_Chimeric.out.junction"),
                    new string[0],
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "fusion_out"));
                return;
            }
            else if (options.Command == "starFusionTest")
            {
                return;
            }

            #endregion STAR Fusion Testing

            #region Proteoform Database Engine

            string proteinDb;
            if (options.SraAccession != null && options.SraAccession.StartsWith("SR"))
            {
                Fastq2ProteinsRunner.RunFromSra(
                    options.BinDirectory,
                    options.BinDirectory,
                    options.Reference,
                    options.Threads,
                    options.SraAccession,
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    options.GeneModelGtfOrGff,
                    out proteinDb);
            }
            else if (options.Fastq1 != null)
            {
                Fastq2ProteinsRunner.RunFromFastqs(
                    options.BinDirectory,
                    options.BinDirectory,
                    options.Reference,
                    options.Threads,
                    options.Fastq2 == null ? new string[] { options.Fastq1 } : new string[] { options.Fastq1, options.Fastq2 },
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    options.GeneModelGtfOrGff,
                    out proteinDb);
            }
            else
            {
                proteinDb = "Error: no fastq or sequence read archive (SRA) was provided.";
            }

            Console.WriteLine("ouptput database to " + proteinDb);
            Console.ReadKey();
            
            #endregion Proteoform Database Engine

        }
    }
}
