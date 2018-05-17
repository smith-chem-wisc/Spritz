using CommandLine;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using ToolWrapperLayer;
using WorkflowLayer;

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

            if (options.Command.Equals(SampleSpecificProteinDBFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                new SnpEffWrapper().DownloadSnpEffDatabase(options.SpritzDirectory, options.AnalysisDirectory, options.Reference);

                if (options.ReferenceVcf == null)
                {
                    var gatk = new GATKWrapper();
                    gatk.DownloadEnsemblKnownVariantSites(options.SpritzDirectory, true, options.Reference);
                    options.ReferenceVcf = gatk.EnsemblKnownSitesPath;
                }

                SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
                flow.Parameters.SpritzDirectory = options.SpritzDirectory;
                flow.Parameters.AnalysisDirectory = options.AnalysisDirectory;
                flow.Parameters.Reference = options.Reference;
                flow.Parameters.Threads = options.Threads;
                flow.Parameters.Fastqs = fastqsSeparated;
                flow.Parameters.StrandSpecific = options.StrandSpecific;
                flow.Parameters.InferStrandSpecificity = options.InferStrandSpecificity;
                flow.Parameters.OverwriteStarAlignment = options.OverwriteStarAlignments;
                flow.Parameters.GenomeStarIndexDirectory = options.GenomeStarIndexDirectory;
                flow.Parameters.GenomeFasta = options.GenomeFasta;
                flow.Parameters.ProteinFasta = options.ProteinFastaPath;
                flow.Parameters.ReferenceGeneModelGtfOrGff = options.GeneModelGtfOrGff;
                flow.Parameters.NewGeneModelGtfOrGff = options.NewGeneModelGtfOrGff;
                flow.Parameters.EnsemblKnownSitesPath = options.ReferenceVcf;
                flow.Parameters.UniProtXmlPath = options.UniProtXml;
                flow.GenerateSAVProteinsFromFastqs();

                Console.WriteLine("output databases to " + String.Join(", and ",
                    flow.VariantAnnotatedProteinXmlDatabases.Concat(flow.VariantAppliedProteinXmlDatabases.Concat(flow.IndelAppliedProteinXmlDatabases))));
            }

            if (options.Command.Equals(LncRNADiscoveryFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                LncRNADiscoveryFlow lnc = new LncRNADiscoveryFlow();
                lnc.Parameters.SpritzDirectory = options.SpritzDirectory;
                lnc.Parameters.AnalysisDirectory = options.AnalysisDirectory;
                lnc.Parameters.Reference = options.Reference;
                lnc.Parameters.Threads = options.Threads;
                lnc.Parameters.Fastqs = fastqsSeparated;
                lnc.Parameters.StrandSpecific = options.StrandSpecific;
                lnc.Parameters.InferStrandSpecificity = options.InferStrandSpecificity;
                lnc.Parameters.OverwriteStarAlignment = options.OverwriteStarAlignments;
                lnc.Parameters.GenomeStarIndexDirectory = options.GenomeStarIndexDirectory;
                lnc.Parameters.GenomeFasta = options.GenomeFasta;
                lnc.Parameters.ProteinFasta = options.ProteinFastaPath;
                lnc.Parameters.GeneModelGtfOrGff = options.GeneModelGtfOrGff;
                lnc.LncRNADiscoveryFromFastqs();
                return;
            }

            if (options.Command.Equals(GeneFusionDiscoveryFlow.Command, StringComparison.InvariantCultureIgnoreCase))
            {
                GeneFusionDiscoveryFlow flow = new GeneFusionDiscoveryFlow();
                flow.Parameters.SpritzDirectory = options.SpritzDirectory;
                flow.Parameters.AnalysisDirectory = options.AnalysisDirectory;
                flow.Parameters.Reference = options.Reference;
                flow.Parameters.Threads = options.Threads;
                flow.Parameters.Fastqs = fastqsSeparated;
                flow.DiscoverGeneFusions();
                return;
            }

            if (options.Command.Equals(TransferModificationsFlow.Command))
            {
                string[] xmls = options.UniProtXml.Split(',');
                TransferModificationsFlow transfer = new TransferModificationsFlow();
                transfer.TransferModifications(options.SpritzDirectory, xmls[0], xmls[1]);
                return;
            }

            if (options.Command.Equals(TranscriptQuantificationFlow.Command, StringComparison.InvariantCultureIgnoreCase))
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