using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class LncRNADiscoveryFlow
        : SpritzFlow
    {
        public const string Command = "lncRNADiscovery";

        public LncRNADiscoveryFlow()
            : base(MyWorkflow.LncRnaDiscovery)
        {
        }

        public LncRNADiscoveryParameters Parameters { get; set; } = new LncRNADiscoveryParameters();
        public string SlnckyOutPrefix { get; private set; }
        public List<string> ReconstructedTranscriptModels { get; private set; } = new List<string>();
        public string MergedTranscriptModel { get; private set; }
        public List<string> IsoformResultPaths { get; private set; } = new List<string>();
        public List<string> GeneResultPaths { get; private set; } = new List<string>();

        /// <summary>
        /// lncRNA discovery from fastq files
        /// </summary>
        public void LncRNADiscoveryFromFastqs()
        {
            // Setup and Alignments
            EnsemblDownloadsWrapper ensemblDownloads = new EnsemblDownloadsWrapper();
            ensemblDownloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            AlignmentFlow alignment = new AlignmentFlow();
            alignment.Parameters = new AlignmentParameters();
            alignment.Parameters.SpritzDirectory = Parameters.SpritzDirectory;
            alignment.Parameters.AnalysisDirectory = Parameters.AnalysisDirectory;
            alignment.Parameters.Reference = Parameters.Reference;
            alignment.Parameters.Threads = Parameters.Threads;
            alignment.Parameters.Fastqs = Parameters.Fastqs;
            alignment.Parameters.ExperimentType = ExperimentType.RNASequencing;
            alignment.Parameters.StrandSpecific = Parameters.StrandSpecific;
            alignment.Parameters.InferStrandSpecificity = Parameters.InferStrandSpecificity;
            alignment.Parameters.OverwriteStarAlignment = Parameters.OverwriteStarAlignment;
            alignment.Parameters.GenomeStarIndexDirectory = Parameters.GenomeStarIndexDirectory;
            alignment.Parameters.ReorderedFastaPath = ensemblDownloads.ReorderedFastaPath;
            alignment.Parameters.GeneModelGtfOrGffPath = Parameters.GeneModelGtfOrGff;
            alignment.Parameters.UseReadSubset = Parameters.UseReadSubset;
            alignment.Parameters.ReadSubset = Parameters.ReadSubset;

            alignment.PerformAlignment();
            ensemblDownloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFasta);
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.AnalysisDirectory, Parameters.GeneModelGtfOrGff, ensemblDownloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.GffOrGtf2Bed12(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, filteredGeneModelForScalpel);

            // Transcript Reconstruction
            StringtieWrapper stringtie = new StringtieWrapper();
            stringtie.TranscriptReconstruction(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, Parameters.GeneModelGtfOrGff, ensemblDownloads.EnsemblGenome,
                Parameters.StrandSpecific, Parameters.InferStrandSpecificity, alignment.SortedBamFiles, true);
            ReconstructedTranscriptModels = stringtie.FilteredTranscriptGtfPaths;
            MergedTranscriptModel = stringtie.FilteredMergedGtfPath;

            // Annotate lncRNAs
            foreach (string gtf in ReconstructedTranscriptModels)
            {
                string slnckyScriptName = WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "SlnckyAnnotation.bash");
                SlnckyOutPrefix = Path.Combine(Path.GetDirectoryName(gtf), Path.GetFileNameWithoutExtension(gtf) + ".slnckyOut", "annotated");
                WrapperUtility.GenerateAndRunScript(slnckyScriptName,
                    SlnckyWrapper.Annotate(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads,
                        gtf, Parameters.Reference, SlnckyOutPrefix)).WaitForExit();
            }

            // Write quantification tables for differential expression analysis (using stringtie TPM values)
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (LncRNADiscoveryParameters)parameters;
            LncRNADiscoveryFromFastqs();
        }
    }
}