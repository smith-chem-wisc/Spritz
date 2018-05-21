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
        public string MergedGtfPath { get; private set; }
        public string MergedFilteredGtfPath { get; private set; }
        public List<string> RsemOutPrefixes { get; private set; } = new List<string>();
        public List<string> ReconstructedTranscriptModels { get; private set; } = new List<string>();
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
                ensemblDownloads.ReorderedFastaPath,
                Parameters.GeneModelGtfOrGff,
                Parameters.UseReadSubset,
                Parameters.ReadSubset);
            alignment.PerformTwoPassAlignment();
            ensemblDownloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFasta);
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.AnalysisDirectory, Parameters.GeneModelGtfOrGff, ensemblDownloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, filteredGeneModelForScalpel);

            // Transcript Reconstruction
            StringTieWrapper stringtie = new StringTieWrapper();
            stringtie.TranscriptReconstruction(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, Parameters.GeneModelGtfOrGff, ensemblDownloads.EnsemblGenome,
                Parameters.StrandSpecific, Parameters.InferStrandSpecificity, alignment.SortedBamFiles);
            ReconstructedTranscriptModels = stringtie.TranscriptGtfPaths;
            MergedGtfPath = stringtie.MergedGtfPath;
            MergedFilteredGtfPath = stringtie.MergedFilteredGtfPath;

            // Transcript Quantification
            foreach (var fastq in Parameters.Fastqs)
            {
                TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
                quantification.Parameters = new TranscriptQuantificationParameters(
                    Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.GenomeFasta, Parameters.Threads, MergedFilteredGtfPath, RSEMAlignerOption.STAR,
                    Parameters.StrandSpecific ? Strandedness.Forward : Strandedness.None,
                    fastq, Parameters.DoOutputQuantificationBam);
                quantification.QuantifyTranscripts();
                RsemOutPrefixes.Add(quantification.RsemOutputPrefix);
                IsoformResultPaths.Add(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix);
                GeneResultPaths.Add(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix);
            }

            // Annotate lncRNAs
            string slnckyScriptName = WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "SlcnkyAnnotation.bash");
            SlnckyOutPrefix = Path.Combine(Path.GetDirectoryName(MergedFilteredGtfPath), Path.GetFileNameWithoutExtension(MergedFilteredGtfPath) + ".slnckyOut", "annotated");
            WrapperUtility.GenerateAndRunScript(slnckyScriptName,
                SlnckyWrapper.Annotate(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads,
                    MergedFilteredGtfPath, Parameters.Reference, SlnckyOutPrefix)).WaitForExit();
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