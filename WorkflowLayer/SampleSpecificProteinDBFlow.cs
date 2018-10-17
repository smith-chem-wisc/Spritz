using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    /// <summary>
    /// Create a database of proteins with single amino acid variants (SAVs) using GATK and SnpEff.
    /// </summary>
    public class SampleSpecificProteinDBFlow
        : SpritzFlow
    {
        public const string Command = "proteins";

        public SampleSpecificProteinDBFlow()
            : base(MyWorkflow.SampleSpecificProteinDB)
        {
        }

        public SampleSpecificProteinDBParameters Parameters { get; set; } = new SampleSpecificProteinDBParameters();
        public List<string> IndelAppliedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> IndelAppliedProteinXmlDatabases { get; private set; } = new List<string>();
        public List<string> VariantAnnotatedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> VariantAnnotatedProteinXmlDatabases { get; private set; } = new List<string>();
        public List<string> VariantAppliedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> VariantAppliedProteinXmlDatabases { get; private set; } = new List<string>();
        private EnsemblDownloadsWrapper Downloads { get; set; } = new EnsemblDownloadsWrapper();

        /// <summary>
        /// Generate sample specific protein database starting with fastq files
        /// </summary>
        public void GenerateSampleSpecificProteinDatabases()
        {
            // Download references and align reads
            Downloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            STARAlignmentFlow alignment = new STARAlignmentFlow();
            if (Parameters.Fastqs != null)
            {
                alignment.Parameters = new STARAlignmentParameters(
                    Parameters.SpritzDirectory,
                    Parameters.AnalysisDirectory,
                    Parameters.Reference,
                    Parameters.Threads,
                    Parameters.Fastqs,
                    Parameters.StrandSpecific,
                    Parameters.InferStrandSpecificity,
                    Parameters.OverwriteStarAlignment,
                    Parameters.GenomeStarIndexDirectory,
                    Downloads.ReorderedFastaPath,
                    Parameters.ReferenceGeneModelGtfOrGff,
                    Parameters.UseReadSubset,
                    Parameters.ReadSubset);
                alignment.PerformTwoPassAlignment();
                Downloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFastaPath);
            }
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.AnalysisDirectory, Parameters.ReferenceGeneModelGtfOrGff, Downloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.GffOrGtf2Bed12(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, filteredGeneModelForScalpel);
            GeneModel referenceGeneModel = new GeneModel(Downloads.EnsemblGenome, Parameters.ReferenceGeneModelGtfOrGff);
            string referenceGeneModelProteinXml = Path.Combine(Path.GetDirectoryName(Parameters.ReferenceGeneModelGtfOrGff), Path.GetFileNameWithoutExtension(Parameters.ReferenceGeneModelGtfOrGff) + ".protein.xml"); // used if no fastqs are provided

            // Merge reference gene model and a new gene model (either specified or stringtie-generated)
            string newGeneModelPath = Parameters.NewGeneModelGtfOrGff;
            string mergedGeneModelWithCdsPath = null;
            string mergedGeneModelProteinXml = null;
            string reference = Parameters.Reference;
            if (Parameters.DoTranscriptIsoformAnalysis)
            {
                StringtieWrapper stringtie = new StringtieWrapper();
                if (newGeneModelPath == null)
                {
                    stringtie.TranscriptReconstruction(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, Parameters.ReferenceGeneModelGtfOrGff, Downloads.EnsemblGenome, Parameters.StrandSpecific, Parameters.InferStrandSpecificity, alignment.SortedBamFiles, true);
                    newGeneModelPath = stringtie.FilteredMergedGtfPath;
                }
                else
                {
                    mergedGeneModelWithCdsPath = Path.Combine(Path.GetDirectoryName(newGeneModelPath), Path.GetFileNameWithoutExtension(newGeneModelPath) + ".withcds.gtf");
                    mergedGeneModelProteinXml = Path.Combine(Path.GetDirectoryName(mergedGeneModelWithCdsPath), Path.GetFileNameWithoutExtension(mergedGeneModelWithCdsPath) + ".protein.xml"); // used if no fastqs are provided, but transcript isoform analysis is performed
                    newGeneModelPath = EnsemblDownloadsWrapper.ConvertFirstColumnUCSC2Ensembl(Parameters.SpritzDirectory, Parameters.Reference, Parameters.NewGeneModelGtfOrGff);
                    string mergedGeneModelPath = Path.Combine(Path.GetDirectoryName(newGeneModelPath), Path.GetFileNameWithoutExtension(newGeneModelPath) + ".merged.gtf");
                    WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "MergeTranscriptModels.bash"),
                        StringtieWrapper.MergeTranscriptPredictions(Parameters.SpritzDirectory, Parameters.ReferenceGeneModelGtfOrGff, new List<string> { newGeneModelPath }, mergedGeneModelPath)).WaitForExit();
                    newGeneModelPath = mergedGeneModelPath;
                }

                GeneModel newGeneModel = new GeneModel(Downloads.EnsemblGenome, newGeneModelPath);

                // Determine CDS from start codons of reference gene model
                // In the future, we could also try ORF finding to expand this (e.g. https://github.com/TransDecoder/TransDecoder/wiki)
                newGeneModel.CreateCDSFromAnnotatedStartCodons(referenceGeneModel);
                newGeneModel.PrintToGTF(mergedGeneModelWithCdsPath);
            }

            // SnpEff databases or outputing protein XMLs from gene models
            if (Parameters.DoTranscriptIsoformAnalysis) // isoform analysis, so generate a new snpeff database
            { 
                reference = SnpEffWrapper.GenerateDatabase(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Downloads.ReorderedFastaPath, Parameters.ProteinFastaPath, mergedGeneModelWithCdsPath);

                if (Parameters.Fastqs == null) // isoform analysis without fastqs, so generate a protein database directly from merged gtf
                    mergedGeneModelProteinXml = SnpEffWrapper.GenerateXmlDatabaseFromReference(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, reference, mergedGeneModelWithCdsPath);
            }
            else if (Parameters.Fastqs != null) // no isoform analysis, but there are are fastqs
            {
                new SnpEffWrapper().DownloadSnpEffDatabase(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Reference);
            }
            else // no isoform analysis and no fastqs
            {
                referenceGeneModelProteinXml = SnpEffWrapper.GenerateXmlDatabaseFromReference(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Reference, Parameters.ReferenceGeneModelGtfOrGff);
            }

            // Variant Calling
            VariantCallingFlow variantCalling = new VariantCallingFlow();
            if (Parameters.Fastqs != null)
            {
                variantCalling.CallVariants(
                    Parameters.SpritzDirectory,
                    Parameters.AnalysisDirectory,
                    reference,
                    Parameters.Threads,
                    sortedBed12Path,
                    Parameters.EnsemblKnownSitesPath,
                    alignment.DedupedBamFiles,
                    Downloads.ReorderedFastaPath,
                    Downloads.EnsemblGenome,
                    Parameters.QuickSnpEffWithoutStats,
                    Parameters.IndelFinder);
            }

            // Gene Fusion Discovery
            List<Protein> fusionProteins = new List<Protein>();
            if (Parameters.DoFusionAnalysis)
            {
                GeneFusionDiscoveryFlow fusion = new GeneFusionDiscoveryFlow();
                fusion.Parameters.SpritzDirectory = Parameters.SpritzDirectory;
                fusion.Parameters.AnalysisDirectory = Parameters.AnalysisDirectory;
                fusion.Parameters.Reference = Parameters.Reference;
                fusion.Parameters.Threads = Parameters.Threads;
                fusion.Parameters.Fastqs = Parameters.Fastqs;
                fusion.DiscoverGeneFusions();
                fusionProteins = fusion.FusionProteins;
            }

            // Transfer features from UniProt
            List<string> xmlsToUse = null;
            if (variantCalling.CombinedAnnotatedProteinXmlPaths.Count > 0)
                xmlsToUse = variantCalling.CombinedAnnotatedProteinXmlPaths;
            else
                xmlsToUse = new List<string> { Parameters.DoTranscriptIsoformAnalysis ? mergedGeneModelProteinXml : referenceGeneModelProteinXml };
            new TransferModificationsFlow().TransferModifications(Parameters.SpritzDirectory, Parameters.UniProtXmlPath, xmlsToUse, fusionProteins);
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (SampleSpecificProteinDBParameters)parameters;
            GenerateSampleSpecificProteinDatabases();
        }
    }
}