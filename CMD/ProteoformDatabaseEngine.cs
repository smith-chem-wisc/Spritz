using System;
using System.Collections.Generic;
using System.Linq;
using WorkflowLayer;
using ToolWrapperLayer;
using CommandLine;
using System.IO;
using System.Reflection;

namespace CMD
{
    class ProteoformDatabaseEngine
    {
        public static void Main(string[] args)
        {

            #region Setup

            if (args.Contains("setup"))
            {
                InstallFlow.Run(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
                return;
            }

            Options options = new Options();
            bool isValid = Parser.Default.ParseArgumentsStrict(args, options);

            #endregion Setup

            #region STAR Fusion Testing

            // TODO make a star fusion 2 protein runner instead of this mess...
            //bool validStarFusionTest = options.AnalysisDirectory != null
            //    && options.BinDirectory != null
            //    && (options.Fastq1 != null || options.st)
            //if (options.Command == "starFusionTest" && !starFusionRequirements.Any(x => x == null))
            //{
            //    Directory.CreateDirectory(Path.Combine(options.AnalysisDirectory, "fusion_out"));
            //    STARFusionWrapper.Install(options.BinDirectory);
            //    STARFusionWrapper.RunStarFusion(options.BinDirectory,
            //        "grch37",
            //        8,
            //        Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR791578_hg19_Chimeric.out.junction"),
            //        new string[0],
            //        Path.Combine(TestContext.CurrentContext.TestDirectory, "fusion_out"));
            //    return;
            //}
            //else if (options.Command == "starFusionTest")
            //{
            //    return;
            //}

            #endregion STAR Fusion Testing

            #region Proteoform Database Engine

            EnsemblDownloadsWrapper.DownloadReferences(
                options.BinDirectory, 
                options.AnalysisDirectory, 
                options.Reference, 
                out string genomeFastaPath,
                out string gtfGeneModelPath,
                out string gff3GeneModelPath);

            if (options.GenomeStarIndexDirectory == null)
                options.GenomeStarIndexDirectory = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath));
            if (options.GenomeFasta == null)
                options.GenomeFasta = genomeFastaPath;
            if (options.GeneModelGtfOrGff == null)
                options.GeneModelGtfOrGff = gff3GeneModelPath;
            if (options.ReferenceVcf == null)
            {
                GATKWrapper.DownloadEnsemblKnownVariantSites(options.BinDirectory, options.AnalysisDirectory, true, options.Reference, out string ensemblVcfPath);
                options.ReferenceVcf = ensemblVcfPath;
            }

            List<string> proteinDatabases;
            if (options.SraAccession != null && options.SraAccession.StartsWith("SR"))
            {
                Fastq2ProteinsEngine.RunFromSra(
                    options.BinDirectory,
                    options.AnalysisDirectory,
                    options.Reference,
                    options.Threads,
                    options.SraAccession,
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    options.GeneModelGtfOrGff,
                    options.ReferenceVcf,
                    out proteinDatabases);
            }
            else if (options.Fastq1 != null)
            {
                Fastq2ProteinsEngine.RunFromFastqs(
                    options.BinDirectory,
                    options.AnalysisDirectory,
                    options.Reference,
                    options.Threads,
                    options.Fastq2 == null ? new string[] { options.Fastq1 } : new string[] { options.Fastq1, options.Fastq2 },
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    options.GeneModelGtfOrGff,
                    options.ReferenceVcf,
                    out proteinDatabases);
            }
            else
            {
                proteinDatabases = new List<string> { "Error: no fastq or sequence read archive (SRA) was provided." };
            }

            Console.WriteLine("output databases to " + String.Join(", and ", proteinDatabases));
            Console.ReadKey();
            
            #endregion Proteoform Database Engine

        }
    }
}
