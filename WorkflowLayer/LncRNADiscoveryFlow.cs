using Proteogenomics;
using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class LncRNADiscoveryFlow
        : SpritzFlow
    {
        public LncRNADiscoveryFlow()
            : base(MyWorkflow.LncRnaDiscovery)
        {
            Parameters = new Parameters();
        }

        public Parameters Parameters { get; set; }
        public string SlnckyOutPrefix { get; private set; }
        public string MergedGtfPath { get; private set; }
        public List<string> RsemOutPrefixes { get; private set; } = new List<string>();
        public List<string> ReconstructedTranscriptModels { get; private set; } = new List<string>();
        public List<string> IsoformResultPaths { get; private set; } = new List<string>();
        public List<string> GeneResultPaths { get; private set; } = new List<string>();

        /// <summary>
        /// Generate sample specific database starting with SRA accession number
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="reference"></param>
        /// <param name="threads"></param>
        /// <param name="sraAccession"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="overwriteStarAlignment"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="proteinFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="ensemblKnownSitesPath"></param>
        /// <param name="proteinVariantDatabases"></param>
        /// <param name="useReadSubset"></param>
        /// <param name="readSubset"></param>
        public void LncRNADiscoveryFromSra(
            string binDirectory, string analysisDirectory, string reference, int threads, string sraAccession,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, bool doOutputQuantificaitonBam,
            bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
                sratoolkit.Fetch(binDirectory, sra, analysisDirectory);
                fastqs.Add(sratoolkit.FastqPaths);
            }
            LncRNADiscoveryFromFastqs(
                binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity,
                overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, doOutputQuantificaitonBam, useReadSubset, readSubset);
        }

        /// <summary>
        /// lncRNA discovery from fastq files
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="reference"></param>
        /// <param name="threads"></param>
        /// <param name="fastqs"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="overwriteStarAlignment"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="proteinFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="doOutputQuantificaitonBam"></param>
        /// <param name="slnckyOutPrefix"></param>
        /// <param name="mergedGtfPath"></param>
        /// <param name="rsemOutPrefixes"></param>
        /// <param name="reconstructedTranscriptModels"></param>
        /// <param name="useReadSubset"></param>
        /// <param name="readSubset"></param>
        public void LncRNADiscoveryFromFastqs(
            string binDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, bool doOutputQuantificaitonBam,
            bool useReadSubset = false, int readSubset = 300000)
        {
            // Setup and Alignments
            STARAlignmentFlow alignment = new STARAlignmentFlow();
            EnsemblDownloadsWrapper.PrepareEnsemblGenomeFasta(genomeFasta, out Genome ensemblGenome, out string reorderedFasta);
            alignment.PerformTwoPassAlignment(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, reorderedFasta, proteinFasta, geneModelGtfOrGff, useReadSubset, readSubset);
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(binDirectory, proteinFasta, out var proteinSequences, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);
            EnsemblDownloadsWrapper.FilterGeneModel(binDirectory, geneModelGtfOrGff, ensemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(binDirectory, filteredGeneModelForScalpel, genomeFasta);

            // Transcript Reconstruction
            StringTieWrapper stringtie = new StringTieWrapper();
            stringtie.TranscriptReconstruction(binDirectory, analysisDirectory, threads, geneModelGtfOrGff, ensemblGenome, strandSpecific, inferStrandSpecificity, alignment.SortedBamFiles);
            ReconstructedTranscriptModels = stringtie.TranscriptGtfPaths;
            MergedGtfPath = stringtie.MergedGtfPath;

            // Transcript Quantification
            foreach (var fastq in fastqs)
            {
                TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
                quantification.QuantifyTranscripts(
                    binDirectory, genomeFasta, threads, MergedGtfPath, RSEMAlignerOption.STAR,
                    strandSpecific ? Strandedness.Forward : Strandedness.None,
                    fastq, doOutputQuantificaitonBam);
                RsemOutPrefixes.Add(quantification.RsemOutputPrefix);
                IsoformResultPaths.Add(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix);
                GeneResultPaths.Add(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix);
            }

            // Annotate lncRNAs
            string slnckyScriptName = Path.Combine(binDirectory, "scripts", "SlcnkyAnnotation.bash");
            SlnckyOutPrefix = Path.Combine(Path.GetDirectoryName(MergedGtfPath), Path.GetFileNameWithoutExtension(MergedGtfPath) + ".slnckyOut", "annotated");
            WrapperUtility.GenerateAndRunScript(slnckyScriptName, SlnckyWrapper.Annotate(binDirectory, analysisDirectory, threads, MergedGtfPath, reference, SlnckyOutPrefix)).WaitForExit();
        }

        public static void Test(string test)
        {
            string script_path = Path.Combine(test, "scripts", "test.sh");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(test),
                "echo " + "HaHa"
            }).WaitForExit();
        }

        protected override void RunSpecific(string OutputFolder, List<string> genomeFastaList, List<string> geneSetList, List<string> rnaSeqFastqList)
        {
            Test(OutputFolder);
        }
    }
}