using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    /// <summary>
    /// Contains a workflow for running a 2-pass STAR alignment on a set of FASTQ files.
    /// </summary>
    public class AlignmentFlow
        : SpritzFlow
    {
        public const string Command = "align";

        public AlignmentFlow()
            : base(MyWorkflow.STARAlignment)
        {
        }

        public AlignmentParameters Parameters { get; set; }
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
        /// Runs a two-pass STAR alignment for a given set of RNA-Seq fastq files,
        /// or it performs a Bowtie2 alignment for WGS or exome sequencing files.
        /// </summary>
        public void PerformAlignment()
        {
            int starThreads = Math.Min(18, Parameters.Threads); // 18 max, otherwise it throws a segmentation fault in sorting the BAM files
            if (Parameters.ExperimentType == ExperimentType.RNASequencing)
            {
                // Alignment preparation
                WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "GenomeGenerate.bash"),
                    STARWrapper.GenerateGenomeIndex(
                        Parameters.SpritzDirectory,
                        Parameters.Threads,
                        Parameters.GenomeStarIndexDirectory,
                        new string[] { Parameters.ReorderedFastaPath },
                        Parameters.GeneModelGtfOrGffPath,
                        Parameters.Fastqs))
                    .WaitForExit();

                // there's trouble with the number of open files for sorting and stuff, which increases with the number of threads
                // 18 is the max that works with the default max number of open files
                TwoPassAlignment(starThreads, Parameters.OverwriteStarAlignment);
            }
            else
            {
                foreach (string[] fastq in Parameters.Fastqs)
                {
                    SkewerWrapper.Trim(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, 19, fastq, false, out string[] trimmedFastqs, out string skewerLog);
                    FastqsForAlignment.Add(trimmedFastqs);
                }
                TopHatWrapper.GenerateBowtieIndex(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.ReorderedFastaPath, out string bowtieIndexPrefix);
                List<string> alignmentCommands = new List<string>();
                foreach (string[] fastq in FastqsForAlignment)
                {
                    // alignment
                    alignmentCommands.AddRange(TopHatWrapper.Bowtie2Align(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, 
                        bowtieIndexPrefix, Parameters.Threads, fastq, Parameters.StrandSpecific, out string sortedBamPath));
                    alignmentCommands.Add(SamtoolsWrapper.IndexBamCommand(sortedBamPath));

                    // mark duplicates
                    GATKWrapper gatk = new GATKWrapper(1);
                    alignmentCommands.AddRange(gatk.PrepareBamAndFasta(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, sortedBamPath, Parameters.ReorderedFastaPath, Parameters.Reference));
                    alignmentCommands.Add(SamtoolsWrapper.IndexBamCommand(gatk.PreparedBamPath));

                    SortedBamFiles.Add(sortedBamPath);
                    DedupedBamFiles.Add(gatk.PreparedBamPath);
                }
                WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "BowtieAlignment.bash"), alignmentCommands).WaitForExit();
            }
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
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "GenomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(spritzDirectory, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff, new List<string[]> { fastqPaths }))
                .WaitForExit();

            STARWrapper.SubsetFastqs(spritzDirectory, analysisDirectory, fastqPaths, 30000, analysisDirectory, out string[] subsetFastqs);

            string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "AlignSubset.bash"),
                STARWrapper.BasicAlignReadCommands(spritzDirectory, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                .WaitForExit();
            BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, new Genome(reorderedFasta), 0.8);
            return bamProperties;
        }

        /// <summary>
        /// Performs the bulk of two-pass alignments
        /// </summary>
        private void TwoPassAlignment(int threads, bool overWriteStarAlignment)
        {
            // Trimming and strand specificity
            Genome genome = new Genome(Parameters.ReorderedFastaPath);
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
                        BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, Parameters.GeneModelGtfOrGffPath, new Genome(Parameters.ReorderedFastaPath), 0.8);
                        localStrandSpecific = bamProperties.Strandedness != Strandedness.None;
                    }
                }

                SkewerWrapper.Trim(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, threads, 19, fqForAlignment, false, out string[] trimmedFastqs, out string skewerLog);
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
            alignmentCommands.AddRange(STARWrapper.GenerateGenomeIndex(Parameters.SpritzDirectory, threads, SecondPassGenomeDirectory, new string[] { Parameters.ReorderedFastaPath }, Parameters.GeneModelGtfOrGffPath, Parameters.Fastqs, spliceJunctionStartDatabase));
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
            Parameters = (AlignmentParameters)parameters;
            PerformAlignment();
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