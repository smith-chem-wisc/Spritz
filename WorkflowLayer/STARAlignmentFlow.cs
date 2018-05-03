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
        : SpritzFlow
    {
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
            WrapperUtility.GenerateAndRunScript(Path.Combine(Parameters.SpritzDirectory, "scripts", "genomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(Parameters.SpritzDirectory, Parameters.Threads, Parameters.GenomeStarIndexDirectory, new string[] { Parameters.ReorderedFasta }, Parameters.GeneModelGtfOrGff))
                .WaitForExit();

            TwoPassAlignment(Parameters.Threads);

            // Rerun if the workaround above changed the number of threads
            int threads = GetDebuggedThreadCount();
            if (threads >= 0 && threads != Parameters.Threads)
            {
                Clear();
                TwoPassAlignment(threads);
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
            WrapperUtility.GenerateAndRunScript(Path.Combine(spritzDirectory, "scripts", "genomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(spritzDirectory, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff))
                .WaitForExit();

            STARWrapper.SubsetFastqs(spritzDirectory, fastqPaths, 30000, analysisDirectory, out string[] subsetFastqs);

            string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
            WrapperUtility.GenerateAndRunScript(Path.Combine(spritzDirectory, "scripts", "alignSubset.bash"),
                STARWrapper.BasicAlignReadCommands(spritzDirectory, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                .WaitForExit();
            BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, new Genome(reorderedFasta), 0.8);
            return bamProperties;
        }

        /// <summary>
        /// STAR has this bug where the Linux environment runs out of room to open more files, see https://github.com/smith-chem-wisc/Spritz/issues/76
        /// It makes these temp files for each thread, numbered starting with 0, and the last one is the one it ran out of room on
        /// </summary>
        /// <returns></returns>
        public int GetDebuggedThreadCount()
        {
            int threads = -1;
            List<string> tmpFilePrefixes = OutputPrefixes.Select(p => Path.Combine(p + "_STARtmp", STARWrapper.ThreadCheckFilePrefix)).ToList();
            foreach (string file in tmpFilePrefixes)
            {
                string tmpdirectory = Path.GetDirectoryName(file);
                if (Directory.Exists(tmpdirectory))
                {
                    var tmpFiles = Directory.GetFiles(tmpdirectory);
                    var tmpFileThreadNumbers = tmpFiles
                        .Where(f => Path.GetFileName(f).StartsWith(STARWrapper.ThreadCheckFilePrefix))
                        .Select(f => Path.GetFileName(f).Substring(STARWrapper.ThreadCheckFilePrefix.Length))
                        .Where(x => int.TryParse(x, out int xx))
                        .Select(x => int.Parse(x))
                        .ToList();
                    threads = tmpFileThreadNumbers.Max() - 1;
                }
            }
            return threads;
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

        /// <summary>
        /// Performs the bulk of two-pass alignments
        /// </summary>
        public void TwoPassAlignment(int threads)
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
                    STARWrapper.SubsetFastqs(Parameters.SpritzDirectory, fqForAlignment, Parameters.ReadSubset, Parameters.AnalysisDirectory, out string[] subsetFastqs);
                    if (Parameters.UseReadSubset)
                    {
                        fqForAlignment = subsetFastqs;
                    }
                    if (Parameters.InferStrandSpecificity)
                    {
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        WrapperUtility.GenerateAndRunScript(Path.Combine(Parameters.SpritzDirectory, "scripts", "alignSubset.bash"),
                            STARWrapper.BasicAlignReadCommands(Parameters.SpritzDirectory, threads, Parameters.GenomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep))
                            .WaitForExit();
                        BAMProperties bamProperties = new BAMProperties(subsetOutPrefix + STARWrapper.BamFileSuffix, Parameters.GeneModelGtfOrGff, new Genome(Parameters.ReorderedFasta), 0.8);
                        localStrandSpecific = bamProperties.Strandedness != Strandedness.None;
                    }
                }

                SkewerWrapper.Trim(Parameters.SpritzDirectory, threads, 19, fqForAlignment, out string[] trimmedFastqs, out string skewerLog);
                fqForAlignment = trimmedFastqs;

                StrandSpecificities.Add(localStrandSpecific);
                FastqsForAlignment.Add(fqForAlignment);
            }

            // Alignment
            List<string> alignmentCommands = new List<string>();
            foreach (string[] fq in FastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SpliceJunctionFileSuffix) || Parameters.OverWriteStarAlignment)
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
                alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(Parameters.SpritzDirectory, threads, SecondPassGenomeDirectory, fq, outPrefix, Parameters.OverWriteStarAlignment, StrandSpecificities[FastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                SortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                DedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
                ChimericSamFiles.Add(outPrefix + STARWrapper.ChimericSamFileSuffix);
                ChimericJunctionFiles.Add(outPrefix + STARWrapper.ChimericJunctionsFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(Parameters.SpritzDirectory, SecondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(Path.Combine(Parameters.SpritzDirectory, "scripts", "alignReads.bash"), alignmentCommands).WaitForExit();
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