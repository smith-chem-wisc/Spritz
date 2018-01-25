using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class STAR2PassAlignFlow
    {
        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string reorderedFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles, out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles, bool useReadSubset = false, int readSubset = 300000)
        {
            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "genomeGenerate.bash"), STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff)).WaitForExit();

            List<string[]> fastqsForAlignment = new List<string[]>();
            List<bool> strandSpecificities = new List<bool>();
            firstPassSpliceJunctions = new List<string>();
            sortedBamFiles = new List<string>();
            dedupedBamFiles = new List<string>();
            chimericSamFiles = new List<string>();
            chimericJunctionFiles = new List<string>();

            // Trimming and strand specificity
            foreach (string[] fq in fastqs)
            {
                SkewerWrapper.Trim(bin, threads, 19, fq, out string[] trimmedFastqs, out string skewerLog);
                string[] fqForAlignment = trimmedFastqs;

                // Infer strand specificity
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
                        WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignSubset.bash"), STARWrapper.BasicAlignReadCommands(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep)).WaitForExit();
                        localStrandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                    }
                }
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
                if (!File.Exists(outPrefix + STARWrapper.SortedBamFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(bin, threads, secondPassGenomeDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));

                }
                sortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                dedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
                chimericSamFiles.Add(outPrefix + STARWrapper.ChimericSamFileSuffix);
                chimericJunctionFiles.Add(outPrefix + STARWrapper.ChimericJunctionsFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, secondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignReads.bash"), alignmentCommands).WaitForExit();
        }
    }
}
