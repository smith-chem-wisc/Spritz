using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace ToolWrapperLayer
{
    public class CufflinksWrapper
    {

        #region Public Properties

        public static string TranscriptsFilename { get; } = "transcripts.gtf";

        public static string SkippedTranscriptsFilename { get; } = "skipped.gtf";

        public static string IsoformAbundanceFilename { get; } = "isoforms.fpkm_tracking";

        public static string GeneAbundanceFilename { get; } = "genes.fpkm_tracking";

        #endregion Public Properties

        #region Public Method

        /// <summary>
        /// Transcript assembly. Note that fragment bias estimation (--frag-bias-correct) and multi-read rescuing (--multi-read-correct) are not used.
        /// These take a lot of time, and they only provide better abundance estimates, which we use RSEM for.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="bamPath"></param>
        /// <param name="outputDirectory"></param>
        /// <param name="geneModelGtfOrGffPath"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        public static void AssembleTranscripts(string binDirectory, int threads, string bamPath, string geneModelGtfOrGffPath, bool strandSpecific, bool inferStrandSpecificity, out string outputDirectory)
        {
            if (inferStrandSpecificity)
            {
                strandSpecific = RSeQCWrapper.CheckStrandSpecificity(binDirectory, bamPath, geneModelGtfOrGffPath, 0.8);
            }

            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksSortCheck");
            outputDirectory = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksOutput");
            string script_name = Path.Combine(binDirectory, "scripts", "cufflinksRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bamPath) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " +
                    "samtools sort -f -@ " + Environment.ProcessorCount.ToString() + " -m " + Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue()) + "M " +
                    WrapperUtility.ConvertWindowsPath(bamPath) + " " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) +
                    "; fi",
                "bam=" +  WrapperUtility.ConvertWindowsPath(bamPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then bam=" + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + "; fi",
                "cufflinks " +
                    " --num-threads " + threads.ToString() +
                    " --GTF-guide " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    (strandSpecific ? "--library-type fr-firststrand" : "") +
                    " $bam",
            }).WaitForExit();
        }

        #endregion Public Method

    }
}
