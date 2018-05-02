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
        public static void LncRNADiscoveryFromSra(
            string binDirectory, string analysisDirectory, string reference, int threads, string sraAccession,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, bool doOutputQuantificaitonBam,
            out string slnckyOutPrefix, out string cuffmergeGtfPath, out List<string> rsemOutPrefixes, out List<string> cufflinksTranscriptModels,
            bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(binDirectory, sra, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            LncRNADiscoveryFromFastqs(
                binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity,
                overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, doOutputQuantificaitonBam,
                out slnckyOutPrefix,
                out cuffmergeGtfPath,
                out rsemOutPrefixes,
                out cufflinksTranscriptModels, useReadSubset, readSubset);
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
        public static void LncRNADiscoveryFromFastqs(
            string binDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, bool doOutputQuantificaitonBam,
            out string slnckyOutPrefix, out string mergedGtfPath, out List<string> rsemOutPrefixes, out List<string> reconstructedTranscriptModels,
            bool useReadSubset = false, int readSubset = 300000)
        {
            // Setup and Alignments
            EnsemblDownloadsWrapper.PrepareEnsemblGenomeFasta(genomeFasta, out Genome ensemblGenome, out string reorderedFasta);
            STARAlignmentFlow.PerformTwoPassAlignment(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, reorderedFasta, proteinFasta, geneModelGtfOrGff, out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles, out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles, useReadSubset, readSubset);
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(binDirectory, proteinFasta, out var proteinSequences, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);
            EnsemblDownloadsWrapper.FilterGeneModel(binDirectory, geneModelGtfOrGff, ensemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(binDirectory, filteredGeneModelForScalpel, genomeFasta);

            // Transcript Reconstruction
            StringTieWrapper.TranscriptReconstruction(binDirectory, analysisDirectory, threads, geneModelGtfOrGff, ensemblGenome, strandSpecific, inferStrandSpecificity,
                sortedBamFiles, out reconstructedTranscriptModels, out mergedGtfPath);

            // Transcript Quantification
            rsemOutPrefixes = new List<string>();
            List<string> isoformResultPaths = new List<string>();
            List<string> geneResultPaths = new List<string>();
            foreach (var fastq in fastqs)
            {
                TranscriptQuantificationFlow.QuantifyTranscripts(
                    binDirectory, genomeFasta, threads, mergedGtfPath, RSEMAlignerOption.STAR,
                    strandSpecific ? Strandedness.Forward : Strandedness.None,
                    fastq, doOutputQuantificaitonBam,
                    out string rsemReferencePrefix, out string rsemOutPrefix);
                rsemOutPrefixes.Add(rsemOutPrefix);
                isoformResultPaths.Add(rsemOutPrefix + RSEMWrapper.IsoformResultsSuffix);
                geneResultPaths.Add(rsemOutPrefix + RSEMWrapper.GeneResultsSuffix);
            }

            // Annotate lncRNAs
            string slnckyScriptName = Path.Combine(binDirectory, "scripts", "SlcnkyAnnotation.bash");
            slnckyOutPrefix = Path.Combine(Path.GetDirectoryName(mergedGtfPath), Path.GetFileNameWithoutExtension(mergedGtfPath) + ".slnckyOut", "annotated");
            WrapperUtility.GenerateAndRunScript(slnckyScriptName, SlnckyWrapper.Annotate(binDirectory, analysisDirectory, threads, mergedGtfPath, reference, slnckyOutPrefix)).WaitForExit();
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