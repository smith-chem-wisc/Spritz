using Proteogenomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using System;

namespace WorkflowLayer
{
    /// <summary>
    /// Contains a workflow for running a 2-pass STAR alignment on a set of FASTQ files.
    /// </summary>
    public class STARAlignmentFlow
        : SpritzFlow
    {
        public const string Command = "a";
    
        public STARAlignmentFlow()
            : base(MyWorkflow.STARAlignment)
        {
        }

        public STARAlignmentParameters Parameters { get; set; }
        public List<string> OutputPrefixes { get; set; } = new List<string>();
        public List<string> FirstPassSpliceJunctions { get; private set; } = new List<string>();
        public string SecondPassGenomeDirectory { get; private set; }
        public List<string> SortedBamFiles { get; private set; } = new List<string>();
        public List<string> DedupedBamFiles { get; private set; } = new List<string>();
        public List<string> ChimericSamFiles { get; private set; } = new List<string>();
        public List<string> ChimericJunctionFiles { get; private set; } = new List<string>();
        public List<string[]> FastqsForAlignment { get; private set; } = new List<string[]>();
        public List<bool> StrandSpecificities { get; private set; } = new List<bool>();

        /// <summary>
        /// Runs a two-pass alignment for a given set of fastq files.
        /// </summary>
        public void PerformTwoPassAlignment()
        {
            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "GenomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(
                    Parameters.SpritzDirectory, 
                    Parameters.Threads, 
                    Parameters.GenomeStarIndexDirectory,
                    new string[] { Parameters.ReorderedFasta }, 
                    Parameters.GeneModelGtfOrGff))
                .WaitForExit();

            // there's trouble with the number of open files for sorting and stuff, which increases with the number of threads
            // 18 is the max that works with the default max number of open files
            TwoPassAlignment(Math.Min(18, Parameters.Threads), Parameters.OverWriteStarAlignment); 
        }

        /// <summary>
        /// Infers the strandedness of reads based on aligning a subset.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="fastqPaths"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="reorderedFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <returns></returns>
        public static BAMProperties InferStrandedness(string spritzDirectory, string analysisDirectory, int threads, string[] fastqPaths, string genomeStarIndexDirectory,
            string reorderedFasta, string geneModelGtfOrGff)
        {
            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "genomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(spritzDirectory, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff))
                .WaitForExit();

            STARWrapper.SubsetFastqs(spritzDirectory, analysisDirectory, fastqPaths, 30000, analysisDirectory, out string[] subsetFastqs);

            string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "alignSubset.bash"),
                STARWrapper.BasicAlignReadCommands(spritzDirectory, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                .WaitForExit();
            BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, new Genome(reorderedFasta), 0.8);
            return bamProperties;
        }

        /// <summary>
        /// Performs the bulk of two-pass alignments
        /// </summary>
        public void TwoPassAlignment(int threads, bool overWriteStarAlignment)
        {
            // Trimming and strand specificity
            Genome genome = new Genome(Parameters.ReorderedFasta);
            foreach (string[] fq in Parameters.Fastqs)
            {
                // Infer strand specificity before trimming because trimming can change read pairings
                string[] fqForAlignment = fq;
                bool localStrandSpecific = Parameters.StrandSpecific;
                if (Parameters.InferStrandSpecificity || Parameters.UseReadSubset)
                {
                    STARWrapper.SubsetFastqs(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, fqForAlignment, 
                        Parameters.ReadSubset, Parameters.AnalysisDirectory, out string[] subsetFastqs);
                    if (Parameters.UseReadSubset)
                    {
                        fqForAlignment = subsetFastqs;
                    }
                    if (Parameters.InferStrandSpecificity)
                    {
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "AlignSubset.bash"),
                            STARWrapper.BasicAlignReadCommands(Parameters.SpritzDirectory, threads, Parameters.GenomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                            .WaitForExit();
                        BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, Parameters.GeneModelGtfOrGff, new Genome(Parameters.ReorderedFasta), 0.8);
                        localStrandSpecific = bamProperties.Strandedness != Strandedness.None;
                    }
                }

                SkewerWrapper.Trim(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, threads, 19, fqForAlignment, out string[] trimmedFastqs, out string skewerLog);
                fqForAlignment = trimmedFastqs;

                StrandSpecificities.Add(localStrandSpecific);
                FastqsForAlignment.Add(fqForAlignment);
            }

            // Alignment
            List<string> alignmentCommands = new List<string>();
            foreach (string[] fq in FastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SpliceJunctionFileSuffix) || overWriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.FirstPassAlignmentCommands(Parameters.SpritzDirectory, threads, Parameters.GenomeStarIndexDirectory, fq, outPrefix, StrandSpecificities[FastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                }
                FirstPassSpliceJunctions.Add(outPrefix + STARWrapper.SpliceJunctionFileSuffix);
            }
            int uniqueSuffix = 1;
            foreach (string f in FastqsForAlignment.SelectMany(f => f))
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(Parameters.SpritzDirectory, Parameters.GenomeStarIndexDirectory));
            alignmentCommands.AddRange(STARWrapper.ProcessFirstPassSpliceCommands(FirstPassSpliceJunctions, uniqueSuffix, out string spliceJunctionStartDatabase));
            SecondPassGenomeDirectory = Parameters.GenomeStarIndexDirectory + "SecondPass" + uniqueSuffix.ToString();
            alignmentCommands.AddRange(STARWrapper.GenerateGenomeIndex(Parameters.SpritzDirectory, threads, SecondPassGenomeDirectory, new string[] { Parameters.ReorderedFasta }, Parameters.GeneModelGtfOrGff, spliceJunctionStartDatabase));
            foreach (string[] fq in FastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                OutputPrefixes.Add(outPrefix);
                alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(Parameters.SpritzDirectory, threads, SecondPassGenomeDirectory, fq, outPrefix, overWriteStarAlignment, StrandSpecificities[FastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                SortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                DedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
                ChimericSamFiles.Add(outPrefix + STARWrapper.ChimericSamFileSuffix);
                ChimericJunctionFiles.Add(outPrefix + STARWrapper.ChimericJunctionsFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(Parameters.SpritzDirectory, SecondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "AlignReads.bash"), alignmentCommands).WaitForExit();
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (STARAlignmentParameters)parameters;
            PerformTwoPassAlignment();
        }

        private void Clear()
        {
            FirstPassSpliceJunctions.Clear();
            SecondPassGenomeDirectory = null;
            SortedBamFiles.Clear();
            DedupedBamFiles.Clear();
            ChimericSamFiles.Clear();
            ChimericJunctionFiles.Clear();
            FastqsForAlignment.Clear();
            StrandSpecificities.Clear();
        }
    }
}