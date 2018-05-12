using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
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

        public SampleSpecificProteinDBParameters Parameters { get; set; }
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
        public void GenerateSAVProteinsFromFastqs()
        {
            // Download references and align reads
            Downloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            STARAlignmentFlow alignment = new STARAlignmentFlow();
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
            Downloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFasta);
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.SpritzDirectory, Parameters.ReferenceGeneModelGtfOrGff, Downloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, filteredGeneModelForScalpel);
            GeneModel referenceGeneModel = new GeneModel(Downloads.EnsemblGenome, Parameters.ReferenceGeneModelGtfOrGff);

            // Merge new and reference gene models, if a new one is specified
            Merge(referenceGeneModel, Parameters.NewGeneModelGtfOrGff);

            // Variant Calling
            VariantCallingFlow variantCalling = new VariantCallingFlow();
            variantCalling.CallVariants(
                Parameters.SpritzDirectory,
                Parameters.AnalysisDirectory,
                Parameters.Reference, 
                Parameters.Threads, 
                sortedBed12Path, 
                Parameters.EnsemblKnownSitesPath, 
                alignment.DedupedBamFiles, 
                Downloads.ReorderedFastaPath);

            // Transfer features from UniProt
            TransferModificationsFlow transfer = new TransferModificationsFlow();
            transfer.TransferModifications(Parameters.SpritzDirectory, Parameters.UniProtXmlPath, variantCalling.CombinedAnnotatedProteinXmlPaths);
        }

        /// <summary>
        /// Read in a new gene model and merge it with this one
        /// </summary>
        /// <param name="alternativeGeneModelPath"></param>
        private void Merge(GeneModel referenceGeneModel, string alternativeGeneModelPath)
        {
            GeneModel newGeneModel = null;
            if (alternativeGeneModelPath != null && File.Exists(alternativeGeneModelPath))
            {
                string newGeneModelPath = EnsemblDownloadsWrapper.ConvertFirstColumnEnsembl2UCSC(Parameters.SpritzDirectory, Parameters.Reference, Parameters.NewGeneModelGtfOrGff);
                newGeneModel = new GeneModel(Downloads.EnsemblGenome, newGeneModelPath);
                newGeneModel.CreateCDSFromAnnotatedStartCodons(referenceGeneModel);
            }
            if (newGeneModel != null)
            {
                referenceGeneModel.Merge(newGeneModel);
            }
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (SampleSpecificProteinDBParameters)parameters;
            GenerateSAVProteinsFromFastqs();
        }
    }
}