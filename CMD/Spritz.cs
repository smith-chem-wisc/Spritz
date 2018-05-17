using CommandLine;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using ToolWrapperLayer;
using WorkflowLayer;
using Proteogenomics;

namespace CMD
{
    public class Spritz
    {
        public static void Main(string[] args)
        {
            if (!WrapperUtility.CheckBashSetup())
            {
                throw new FileNotFoundException("The Windows Subsystem for Windows has not been enabled. Please see https://smith-chem-wisc.github.io/Spritz/ for more details.");
            }

            // main setup involves installing tools
            if (args.Contains(ManageToolsFlow.Command))
            {
                ManageToolsFlow.Install(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
                return;
            }

            Parsed<Options> result = Parser.Default.ParseArguments<Options>(args) as Parsed<Options>;
            Options options = result.Value;

            FinishSetup(options);

            bool useSraMethod = options.SraAccession != null && options.SraAccession.StartsWith("SR");
            List<string[]> fastqsSeparated = useSraMethod ?
                SRAToolkitWrapper.GetFastqsFromSras(options.SpritzDirectory, options.AnalysisDirectory, options.SraAccession) :
                SeparateFastqs(options.Fastq1, options.Fastq2);

            #region STAR Fusion Testing

            // TODO make a star fusion 2 protein runner instead of this mess...
            //bool validStarFusionTest = options.AnalysisDirectory != null
            //    && options.SpritzDirectory != null
            //    && (options.Fastq1 != null || options.st)
            //if (options.Command == "starFusionTest" && !starFusionRequirements.Any(x => x == null))
            //{
            //    Directory.CreateDirectory(Path.Combine(options.AnalysisDirectory, "fusion_out"));
            //    STARFusionWrapper.Install(options.SpritzDirectory);
            //    STARFusionWrapper.RunStarFusion(options.SpritzDirectory,
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

            if (options.Command.Equals(LncRNADiscoveryFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                LncRNADiscoveryFlow lncRNAdiscovery = new LncRNADiscoveryFlow();
                lncRNAdiscovery.Parameters = new LncRNADiscoveryParameters(
                        options.SpritzDirectory,
                        options.AnalysisDirectory,
                        options.Reference,
                        options.Threads,
                        fastqsSeparated,
                        options.StrandSpecific,
                        options.InferStrandSpecificity,
                        options.OverwriteStarAlignments,
                        options.GenomeStarIndexDirectory,
                        options.GenomeFasta,
                        options.ProteinFastaPath,
                        options.GeneModelGtfOrGff,
                        true);
                lncRNAdiscovery.LncRNADiscoveryFromFastqs();
                return;
            }

            #endregion lncRNA Discovery Workflow

            #region Infering Strandedness

            if (options.Command.Equals("strandedness"))
            {
                string[] fastqs = options.Fastq2 == null ?
                    new[] { options.Fastq1 } :
                    new[] { options.Fastq1, options.Fastq2 };
                BAMProperties b = STARAlignmentFlow.InferStrandedness(options.SpritzDirectory, options.AnalysisDirectory, options.Threads,
                        fastqs, options.GenomeStarIndexDirectory, options.GenomeFasta, options.GeneModelGtfOrGff);
                Console.WriteLine(b.ToString());
                return;
            }

            #endregion Infering Strandedness

            #region Infering Strandedness

            if (options.Command.Equals(TransferModificationsFlow.Command))
            {
                string[] xmls = options.UniProtXml.Split(',');
                TransferModificationsFlow transfer = new TransferModificationsFlow();
                transfer.TransferModifications(options.SpritzDirectory, xmls[0], xmls[1]);
                return;
            }

            #endregion Infering Strandedness

            #region Transcript Quantification

            if (options.Command.Equals(LncRNADiscoveryFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                foreach (string[] fastq in fastqsSeparated)
                {
                    Strandedness strandedness = options.StrandSpecific ? Strandedness.Forward : Strandedness.None;
                    if (options.InferStrandSpecificity)
                    {
                        var bamProps = STARAlignmentFlow.InferStrandedness(options.SpritzDirectory, options.AnalysisDirectory, options.Threads,
                            fastq, options.GenomeStarIndexDirectory, options.GenomeFasta, options.GeneModelGtfOrGff);
                        strandedness = bamProps.Strandedness;
                    }
                    TranscriptQuantificationFlow quantify = new TranscriptQuantificationFlow();
                    quantify.Parameters = new TranscriptQuantificationParameters(
                        options.SpritzDirectory,
                        options.AnalysisDirectory,
                        options.GenomeFasta,
                        options.Threads,
                        options.GeneModelGtfOrGff,
                        RSEMAlignerOption.STAR,
                        strandedness,
                        fastq,
                        true);
                    quantify.QuantifyTranscripts();
                }
                return;
            }

            #endregion Transcript Quantification

            #region Proteoform Database Engine

            if (options.Command.Equals(SampleSpecificProteinDBFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                new SnpEffWrapper().DownloadSnpEffDatabase(options.SpritzDirectory, options.AnalysisDirectory, options.Reference);

                if (options.ReferenceVcf == null)
                {
                    var gatk = new GATKWrapper();
                    gatk.DownloadEnsemblKnownVariantSites(options.SpritzDirectory, true, options.Reference);
                    options.ReferenceVcf = gatk.EnsemblKnownSitesPath;
                }

                // run the program
                SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
                ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                    options.SpritzDirectory,
                    options.AnalysisDirectory,
                    options.Reference,
                    options.Threads,
                    fastqsSeparated,
                    options.StrandSpecific,
                    options.InferStrandSpecificity,
                    options.OverwriteStarAlignments,
                    options.GenomeStarIndexDirectory,
                    options.GenomeFasta,
                    options.ProteinFastaPath,
                    options.GeneModelGtfOrGff,
                    options.ReferenceVcf,
                    options.UniProtXml);

                ssdbf.GenerateSAVProteinsFromFastqs();

                Console.WriteLine("output databases to " + String.Join(", and ",
                    ssdbf.VariantAnnotatedProteinXmlDatabases.Concat(ssdbf.VariantAppliedProteinXmlDatabases.Concat(ssdbf.IndelAppliedProteinXmlDatabases))));
            }

            #endregion Proteoform Database Engine
        }

        /// <summary>
        /// Always download reference that aren't present and set default options
        /// </summary>
        /// <param name="options"></param>
        public static void FinishSetup(Options options)
        {
            // Download ensembl references and set default paths
            EnsemblDownloadsWrapper downloadsWrapper = new EnsemblDownloadsWrapper();
            downloadsWrapper.DownloadReferences(options.SpritzDirectory, options.SpritzDirectory, options.Reference);
            options.GenomeFasta = options.GenomeFasta ?? downloadsWrapper.GenomeFastaPath;
            options.GeneModelGtfOrGff = options.GeneModelGtfOrGff ?? downloadsWrapper.Gff3GeneModelPath;
            options.GenomeStarIndexDirectory = options.GenomeStarIndexDirectory ?? STARWrapper.GetGenomeStarIndexDirectoryPath(options.GenomeFasta, options.GeneModelGtfOrGff);
            options.ProteinFastaPath = options.ProteinFastaPath ?? downloadsWrapper.ProteinFastaPath;
        }

        /// <summary>
        /// Split the fastq lists into a list of paired strings
        /// </summary>
        /// <param name="fastq1string"></param>
        /// <param name="fastq2string"></param>
        /// <returns></returns>
        private static List<string[]> SeparateFastqs(string fastq1string, string fastq2string)
        {
            List<string[]> fastqsSeparated = null;
            if (fastq1string != null)
            {
                // Parse comma-separated fastq lists
                if (fastq2string != null && fastq1string.Count(x => x == ',') != fastq1string.Count(x => x == ','))
                {
                    throw new ArgumentException("Error: There are a different number of first-strand and second-strand fastq files.");
                }
                string[] fastqs1 = fastq1string.Split(',');
                fastqsSeparated = fastq2string == null ?
                    fastqs1.Select(x => new string[] { x }).ToList() :
                    fastqs1.Select(x => new string[] { x, fastq2string.Split(',')[fastqs1.ToList().IndexOf(x)] }).ToList();
            }
            return fastqsSeparated;
        }
    }
}