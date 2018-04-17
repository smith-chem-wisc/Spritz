using Proteogenomics;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Cufflinks is a program that performs transcript reconstruction from RNA-Seq data. It's particularly useful for genome-guided reconstruction.
    /// </summary>
    public class CufflinksWrapper :
        IInstallable
    {
        #region Public Properties

        /// <summary>
        /// Output filename for reconstructed transcripts gene model.
        /// </summary>
        public static string TranscriptsFilename { get; } = "transcripts.gtf";

        /// <summary>
        /// Output file for skipped transcripts.
        /// </summary>
        public static string SkippedTranscriptsFilename { get; } = "skipped.gtf";

        /// <summary>
        /// Output file for isoform abundance (although I favor using RSEM).
        /// </summary>
        public static string IsoformAbundanceFilename { get; } = "isoforms.fpkm_tracking";

        /// <summary>
        /// Output file for gene abundance (although I favor using RSEM).
        /// </summary>
        public static string GeneAbundanceFilename { get; } = "genes.fpkm_tracking";

        #endregion Public Properties

        #region Installation Methods

        /// <summary>
        /// Writes a script for installing cufflinks.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installCufflinks.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d cufflinks-2.2.1 ]; then",
                "  wget --no-check http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  tar -xvf cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  rm cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  mv cufflinks-2.2.1.Linux_x86_64 cufflinks-2.2.1",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing cufflinks.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Method

        /// <summary>
        /// Transcript assembly. Note that fragment bias estimation (--frag-bias-correct) and multi-read rescuing (--multi-read-correct) are not used.
        /// These take a lot of time, and they only provide better abundance estimates, which we use RSEM for.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="bamPath"></param>
        /// <param name="geneModelGtfOrGffPath"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="outputDirectory"></param>
        public static void AssembleTranscripts(string binDirectory, int threads, string bamPath, string geneModelGtfOrGffPath, Genome genome, bool strandSpecific, bool inferStrandSpecificity, out string outputDirectory)
        {
            bool isStranded = strandSpecific;
            if (inferStrandSpecificity)
            {
                BAMProperties bamProperties = new BAMProperties(bamPath, geneModelGtfOrGffPath, genome, 0.8);
                isStranded = bamProperties.Strandedness != Strandedness.None;
            }

            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksSortCheck");
            outputDirectory = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksOutput");
            string script_name = Path.Combine(binDirectory, "scripts", "cufflinksRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bamPath) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " + SamtoolsWrapper.SortBam(binDirectory, bamPath) + "; fi",
                "bam=" +  WrapperUtility.ConvertWindowsPath(bamPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then bam=" + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + "; fi",
                "cufflinks-2.2.1/cufflinks " +
                    " --num-threads " + threads.ToString() +
                    " --GTF-guide " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    (isStranded ? "--library-type fr-firststrand" : "") +
                    " $bam",
            }).WaitForExit();
        }

        #endregion Public Method
    }
}