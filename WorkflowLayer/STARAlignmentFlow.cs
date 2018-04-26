using Proteogenomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    /// <summary>
    /// Contains a workflow for running a 2-pass STAR alignment on a set of FASTQ files.
    /// </summary>
    public class STARAlignmentFlow
    {
        /// <summary>
        /// Runs a two-pass alignment for a given set of fastq files.
        /// </summary>
        /// <param name="bin"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="reference"></param>
        /// <param name="threads"></param>
        /// <param name="fastqs">List of fastq files, which may be paired in a string[].</param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="overwriteStarAlignment"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="reorderedFasta"></param>
        /// <param name="proteinFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="firstPassSpliceJunctions"></param>
        /// <param name="secondPassGenomeDirectory"></param>
        /// <param name="sortedBamFiles"></param>
        /// <param name="dedupedBamFiles"></param>
        /// <param name="chimericSamFiles"></param>
        /// <param name="chimericJunctionFiles"></param>
        /// <param name="useReadSubset"></param>
        /// <param name="readSubset"></param>
        public static void PerformTwoPassAlignment(string bin, string analysisDirectory, string reference, int threads,
            List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string reorderedFasta, string proteinFasta, string geneModelGtfOrGff,
            out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles,
            out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles,
            bool useReadSubset = false, int readSubset = 30000)
        {
            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "genomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff))
                .WaitForExit();

            List<string[]> fastqsForAlignment = new List<string[]>();
            List<bool> strandSpecificities = new List<bool>();
            firstPassSpliceJunctions = new List<string>();
            sortedBamFiles = new List<string>();
            dedupedBamFiles = new List<string>();
            chimericSamFiles = new List<string>();
            chimericJunctionFiles = new List<string>();

            // Trimming and strand specificity
            Genome genome = new Genome(reorderedFasta);
            foreach (string[] fq in fastqs)
            {
                // Infer strand specificity before trimming because trimming can change read pairings
                string[] fqForAlignment = fq;
                bool localStrandSpecific = strandSpecific;
                if (inferStrandSpecificity || useReadSubset)
                {
                    STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                    if (useReadSubset)
                    {
                        fqForAlignment = subsetFastqs;
                    }
                    if (inferStrandSpecificity)
                    {
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignSubset.bash"),
                            STARWrapper.BasicAlignReadCommands(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                            .WaitForExit();
                        BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, new Genome(reorderedFasta), 0.8);
                        localStrandSpecific = bamProperties.Strandedness != Strandedness.None;
                    }
                }

                SkewerWrapper.Trim(bin, threads, 19, fqForAlignment, out string[] trimmedFastqs, out string skewerLog);
                fqForAlignment = trimmedFastqs;

                strandSpecificities.Add(localStrandSpecific);
                fastqsForAlignment.Add(fqForAlignment);
            }

            // Alignment
            List<string> alignmentCommands = new List<string>();
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SpliceJunctionFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.FirstPassAlignmentCommands(bin, threads, genomeStarIndexDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                }
                firstPassSpliceJunctions.Add(outPrefix + STARWrapper.SpliceJunctionFileSuffix);
            }
            int uniqueSuffix = 1;
            foreach (string f in fastqsForAlignment.SelectMany(f => f))
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, genomeStarIndexDirectory));
            alignmentCommands.AddRange(STARWrapper.ProcessFirstPassSpliceCommands(firstPassSpliceJunctions, uniqueSuffix, out string spliceJunctionStartDatabase));
            secondPassGenomeDirectory = genomeStarIndexDirectory + "SecondPass" + uniqueSuffix.ToString();
            alignmentCommands.AddRange(STARWrapper.GenerateGenomeIndex(bin, threads, secondPassGenomeDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff, spliceJunctionStartDatabase));
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(bin, threads, secondPassGenomeDirectory, fq, outPrefix, overwriteStarAlignment, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                sortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                dedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
                chimericSamFiles.Add(outPrefix + STARWrapper.ChimericSamFileSuffix);
                chimericJunctionFiles.Add(outPrefix + STARWrapper.ChimericJunctionsFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, secondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignReads.bash"), alignmentCommands).WaitForExit();
        }

        /// <summary>
        /// Infers the strandedness of reads based on aligning a subset.
        /// </summary>
        /// <param name="bin"></param>
        /// <param name="threads"></param>
        /// <param name="fastqPaths"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="reorderedFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <returns></returns>
        public static BAMProperties InferStrandedness(string binDirectory, string analysisDirectory, int threads, string[] fastqPaths, string genomeStarIndexDirectory,
            string reorderedFasta, string geneModelGtfOrGff)
        {
            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(Path.Combine(binDirectory, "scripts", "genomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(binDirectory, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff))
                .WaitForExit();

            STARWrapper.SubsetFastqs(binDirectory, fastqPaths, 30000, analysisDirectory, out string[] subsetFastqs);

            string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
            WrapperUtility.GenerateAndRunScript(Path.Combine(binDirectory, "scripts", "alignSubset.bash"),
                STARWrapper.BasicAlignReadCommands(binDirectory, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                .WaitForExit();
            BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, new Genome(reorderedFasta), 0.8);
            return bamProperties;
        }
    }
}