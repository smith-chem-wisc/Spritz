using CommandLine;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using ToolWrapperLayer;
using WorkflowLayer;

namespace CMD
{
    internal class Spritz
    {
        public static void Main(string[] args)
        {
            #region Setup

            // main setup involves installing tools
            if (args.Contains("setup"))
            {
                ManageToolsFlow.Install(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
                return;
            }

            Parsed<Options> result = Parser.Default.ParseArguments<Options>(args) as Parsed<Options>;
            Options options = result.Value;

            // always download references that aren't present
            EnsemblDownloadsWrapper.DownloadReferences(
                options.BinDirectory,
                options.BinDirectory,
                options.Reference,
                out string genomeFastaPath,
                out string gtfGeneModelPath,
                out string gff3GeneModelPath,
                out string proteinFastaPath);

            if (options.GenomeStarIndexDirectory == null)
            {
                options.GenomeStarIndexDirectory = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath));
            }

            if (options.GenomeFasta == null)
            {
                options.GenomeFasta = genomeFastaPath;
            }

            if (options.GeneModelGtfOrGff == null)
            {
                options.GeneModelGtfOrGff = gff3GeneModelPath;
            }

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

            #region lncRNA Discovery Workflow

            if (options.Command.Equals("lncRNADiscovery", StringComparison.InvariantCultureIgnoreCase))
            {
                if (options.Fastq2 != null && options.Fastq1.Count(x => x == ',') != options.Fastq2.Count(x => x == ','))
                {
                    throw new ArgumentException("Unequal count of fastq1 and fastq2 files.");
                }

                if (options.SraAccession != null && options.SraAccession.StartsWith("SR"))
                {
                    LncRNADiscoveryFlow.LncRNADiscoveryFromSra(
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
                        proteinFastaPath,
                        options.GeneModelGtfOrGff,
                        true,
                        out string slnckyOutPrefix,
                        out string cuffmergeGtfPath,
                        out List<string> rsemOutPrefixes,
                        out List<string> cufflinksTranscriptModels);
                }
                else if (options.Fastq1 != null)
                {
                    string[] fastqs1 = options.Fastq1.Split(',');
                    List<string[]> fastqsSeparated = options.Fastq2 == null ?
                        fastqs1.Select(x => new string[] { x }).ToList() :
                        fastqs1.Select(x => new string[] { x, options.Fastq2.Split(',')[fastqs1.ToList().IndexOf(x)] }).ToList();

                    LncRNADiscoveryFlow.LncRNADiscoveryFromFastqs(
                         options.BinDirectory,
                         options.AnalysisDirectory,
                         options.Reference,
                         options.Threads,
                         fastqsSeparated,
                         options.StrandSpecific,
                         options.InferStrandSpecificity,
                         options.OverwriteStarAlignments,
                         options.GenomeStarIndexDirectory,
                         options.GenomeFasta,
                         proteinFastaPath,
                         options.GeneModelGtfOrGff,
                         true,
                         out string slnckyOutPrefix,
                         out string cuffmergeGtfPath,
                         out List<string> rsemOutPrefixes,
                         out List<string> cufflinksTranscriptModels);
                }

                return;
            }

            #endregion lncRNA Discovery Workflow

            #region Infering Strandedness

            if (options.Command.Equals("strandedness"))
            {
                if (options.Fastq2 != null && options.Fastq1.Count(x => x == ',') != options.Fastq2.Count(x => x == ','))
                {
                    throw new ArgumentException("Unequal count of fastq1 and fastq2 files.");
                }

                string[] fastqs = options.Fastq2 == null ? new[] { options.Fastq1 } : new[] { options.Fastq1, options.Fastq2 };
                BAMProperties b = STARAlignmentFlow.InferStrandedness(options.BinDirectory, options.AnalysisDirectory, options.Threads,
                        fastqs, options.GenomeStarIndexDirectory, options.GenomeFasta, options.GeneModelGtfOrGff);
                Console.WriteLine(b.ToString());
                return;
            }

            #endregion Infering Strandedness

            #region Transcript Quantification

            if (options.Command.Equals("quantify", StringComparison.InvariantCultureIgnoreCase))
            {
                if (options.Fastq2 != null && options.Fastq1.Count(x => x == ',') != options.Fastq2.Count(x => x == ','))
                {
                    throw new ArgumentException("Unequal count of fastq1 and fastq2 files.");
                }

                Strandedness strandedness = options.StrandSpecific ? Strandedness.Forward : Strandedness.None;

                if (options.SraAccession != null && options.SraAccession.StartsWith("SR"))
                {
                    TranscriptQuantificationFlow.QuantifyTranscriptsFromSra(
                        options.BinDirectory,
                        options.AnalysisDirectory,
                        options.GenomeFasta,
                        options.Threads,
                        options.GeneModelGtfOrGff,
                        RSEMAlignerOption.STAR,
                        strandedness,
                        options.SraAccession,
                        true,
                        out string referencePrefix,
                        out string outputPrefix);
                }
                else if (options.Fastq1 != null)
                {
                    string[] fastqs = options.Fastq2 == null ? new[] { options.Fastq1 } : new[] { options.Fastq1, options.Fastq2 };
                    if (options.InferStrandSpecificity)
                    {
                        var bamProps = STARAlignmentFlow.InferStrandedness(options.BinDirectory, options.AnalysisDirectory, options.Threads,
                            fastqs, options.GenomeStarIndexDirectory, options.GenomeFasta, options.GeneModelGtfOrGff);
                        strandedness = bamProps.Strandedness;
                    }

                    TranscriptQuantificationFlow.QuantifyTranscripts(
                        options.BinDirectory,
                        options.GenomeFasta,
                        options.Threads,
                        options.GeneModelGtfOrGff,
                        RSEMAlignerOption.STAR,
                        strandedness,
                        fastqs,
                        true,
                        out string referencePrefix,
                        out string outputPrefix);
                }

                return;
            }

            #endregion Transcript Quantification

            #region Proteoform Database Engine

            SnpEffWrapper.DownloadSnpEffDatabase(
                options.BinDirectory,
                options.Reference,
                out string snpEffDatabaseListPath);

            if (options.ReferenceVcf == null)
            {
                GATKWrapper.DownloadEnsemblKnownVariantSites(options.BinDirectory, options.AnalysisDirectory, true, options.Reference, out string ensemblVcfPath);
                options.ReferenceVcf = ensemblVcfPath;
            }

            // run the program

            List<string> proteinDatabases = new List<string>();

            if (options.SraAccession != null && options.SraAccession.StartsWith("SR"))
            {
                SAVProteinDBFlow.GenerateSAVProteinsFromSra(
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
                    proteinFastaPath,
                    options.GeneModelGtfOrGff,
                    options.ReferenceVcf,
                    out proteinDatabases);
            }
            else if (options.Fastq1 != null)
            {
                // Parse comma-separated fastq lists
                if (options.Fastq2 != null && options.Fastq1.Count(x => x == ',') != options.Fastq2.Count(x => x == ','))
                    return;

                string[] fastqs1 = options.Fastq1.Split(',');
                List<string[]> fastqsSeparated = options.Fastq2 == null ?
                    fastqs1.Select(x => new string[] { x }).ToList() :
                    fastqs1.Select(x => new string[] { x, options.Fastq2.Split(',')[fastqs1.ToList().IndexOf(x)] }).ToList();

                SAVProteinDBFlow.GenerateSAVProteinsFromFastqs(
                    options.BinDirectory,
                    options.AnalysisDirectory,
                    options.Reference,
                    options.Threads,
                    fastqsSeparated,
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    proteinFastaPath,
                    options.GeneModelGtfOrGff,
                    options.ReferenceVcf,
                    out proteinDatabases);
            }
            else if (args.Contains("vcf2protein"))
            {
                Genome genome = new Genome(options.GenomeFasta);
                EnsemblDownloadsWrapper.GetImportantProteinAccessions(options.BinDirectory, proteinFastaPath, out var proteinSequences, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);
                GeneModel geneModel = new GeneModel(genome, options.GeneModelGtfOrGff);
                proteinDatabases.Add(SAVProteinDBFlow.WriteSampleSpecificFasta(options.ReferenceVcf, genome, geneModel, options.Reference, proteinSequences, badProteinAccessions, selenocysteineContainingAccessions, 7, Path.Combine(Path.GetDirectoryName(options.ReferenceVcf), Path.GetFileNameWithoutExtension(options.ReferenceVcf))));
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