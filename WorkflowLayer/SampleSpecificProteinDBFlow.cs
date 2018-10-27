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
        public AlignmentFlow Alignment { get; } = new AlignmentFlow();
        public VariantCallingFlow VariantCalling { get; } = new VariantCallingFlow();
        GeneFusionDiscoveryFlow Fusion = new GeneFusionDiscoveryFlow();
        public List<string> VariantAnnotatedProteinXmlDatabases { get; private set; } = new List<string>();
        private EnsemblDownloadsWrapper Downloads { get; set; } = new EnsemblDownloadsWrapper();

        /// <summary>
        /// Generate sample specific protein database starting with fastq files
        /// </summary>
        public void GenerateSampleSpecificProteinDatabases()
        {
            // Download references and align reads
            Downloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            if (Parameters.Fastqs != null && Parameters.ExperimentType.Equals(ExperimentType.RNASequencing))
            {
                Alignment.Parameters = new AlignmentParameters();
                Alignment.Parameters.SpritzDirectory = Parameters.SpritzDirectory;
                Alignment.Parameters.AnalysisDirectory = Parameters.AnalysisDirectory;
                Alignment.Parameters.Reference = Parameters.Reference;
                Alignment.Parameters.Threads = Parameters.Threads;
                Alignment.Parameters.Fastqs = Parameters.Fastqs;
                Alignment.Parameters.ExperimentType = Parameters.ExperimentType;
                Alignment.Parameters.StrandSpecific = Parameters.StrandSpecific;
                Alignment.Parameters.InferStrandSpecificity = Parameters.InferStrandSpecificity;
                Alignment.Parameters.OverwriteStarAlignment = Parameters.OverwriteStarAlignment;
                Alignment.Parameters.GenomeStarIndexDirectory = Parameters.GenomeStarIndexDirectory;
                Alignment.Parameters.ReorderedFastaPath = Downloads.ReorderedFastaPath;
                Alignment.Parameters.GeneModelGtfOrGffPath = Parameters.ReferenceGeneModelGtfOrGff;
                Alignment.Parameters.UseReadSubset = Parameters.UseReadSubset;
                Alignment.Parameters.ReadSubset = Parameters.ReadSubset;

                Alignment.PerformAlignment();
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
                    stringtie.TranscriptReconstruction(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, Parameters.ReferenceGeneModelGtfOrGff, Downloads.EnsemblGenome, Parameters.StrandSpecific, Parameters.InferStrandSpecificity, Alignment.SortedBamFiles, true);
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

            // Gene Fusion Discovery
            List<Protein> fusionProteins = new List<Protein>();
            if (Parameters.DoFusionAnalysis)
            {
                Fusion.Parameters.SpritzDirectory = Parameters.SpritzDirectory;
                Fusion.Parameters.AnalysisDirectory = Parameters.AnalysisDirectory;
                Fusion.Parameters.Reference = Parameters.Reference;
                Fusion.Parameters.Threads = Parameters.Threads;
                Fusion.Parameters.Fastqs = Parameters.Fastqs;
                Fusion.DiscoverGeneFusions();
                fusionProteins = Fusion.FusionProteins;
            }

            // Variant Calling
            if (Parameters.Fastqs != null && !Parameters.SkipVariantAnalysis)
            {
                VariantCalling.CallVariants(
                    Parameters.SpritzDirectory,
                    Parameters.AnalysisDirectory,
                    Parameters.ExperimentType,
                    reference,
                    Parameters.Threads,
                    sortedBed12Path,
                    Parameters.EnsemblKnownSitesPath,
                    Alignment.DedupedBamFiles,
                    Downloads.ReorderedFastaPath,
                    Downloads.EnsemblGenome,
                    Parameters.QuickSnpEffWithoutStats,
                    Parameters.IndelFinder);
            }

            // Transfer features from UniProt
            List<string> xmlsToUse = null;
            if (VariantCalling.CombinedAnnotatedProteinXmlPaths.Count > 0)
                xmlsToUse = VariantCalling.CombinedAnnotatedProteinXmlPaths; 
            // keep, since it might be useful for making a final database: .Concat(new[] { VariantCalling.CombinedAnnotatedProteinXmlPath }).ToList()
            else
                xmlsToUse = new List<string> { Parameters.DoTranscriptIsoformAnalysis ? mergedGeneModelProteinXml : referenceGeneModelProteinXml };
            VariantAnnotatedProteinXmlDatabases = new TransferModificationsFlow().TransferModifications(Parameters.SpritzDirectory, Parameters.UniProtXmlPath, xmlsToUse, fusionProteins);
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