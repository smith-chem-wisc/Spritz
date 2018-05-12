using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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

        /// <summary>
        /// Output filename for cufflinks
        /// </summary>
        public static string CufflinksMergedFilename { get; } = "merged.gtf";

        #endregion Public Properties

        #region Installation Methods

        /// <summary>
        /// Writes a script for installing cufflinks.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallCufflinks.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d cufflinks-2.2.1 ]; then",
                "  wget --no-check http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  tar -xvf cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  rm cufflinks-2.2.1.Linux_x86_64.tar.gz",
                "  mv cufflinks-2.2.1.Linux_x86_64 cufflinks-2.2.1",
                "  cp cufflinks-2.2.1/* /usr/local/bin",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing cufflinks.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Method

        /// <summary>
        /// Transcript assembly. Note that fragment bias estimation (--frag-bias-correct) and multi-read rescuing (--multi-read-correct) are not used.
        /// These take a lot of time, and they only provide better abundance estimates, which we use RSEM for.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="bamPath"></param>
        /// <param name="geneModelGtfOrGffPath"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="outputDirectory"></param>
        public static List<string> AssembleTranscripts(string spritzDirectory, string analysisDirectory, int threads, string bamPath, string geneModelGtfOrGffPath, Genome genome, bool strandSpecific, bool inferStrandSpecificity, out string outputDirectory)
        {
            bool isStranded = strandSpecific;
            if (inferStrandSpecificity)
            {
                BAMProperties bamProperties = new BAMProperties(bamPath, geneModelGtfOrGffPath, genome, 0.8);
                isStranded = bamProperties.Strandedness != Strandedness.None;
            }

            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksSortCheck");
            outputDirectory = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksOutput");
            string script_name = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "CufflinksRun.bash");
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bamPath) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " + SamtoolsWrapper.SortBam(spritzDirectory, bamPath) + "; fi",
                "bam=" +  WrapperUtility.ConvertWindowsPath(bamPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then bam=" + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + "; fi",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(outputDirectory, TranscriptsFilename)) + " || ! -s " + WrapperUtility.ConvertWindowsPath(Path.Combine(outputDirectory, TranscriptsFilename)) + " ]]; then " +
                    "cufflinks-2.2.1/cufflinks " +
                    " --num-threads " + threads.ToString() +
                    " --GTF-guide " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    (isStranded ? "--library-type fr-firststrand" : "") +
                    " $bam" +
                "; fi",
            };
        }

        /// <summary>
        /// Removes transcripts with zero abundance predictions
        /// </summary>
        /// <returns></returns>
        public static List<string> RemoveZeroAbundanceCufflinksPredictionsCommand(string spritzDirectory, string transcriptGtfPath, out string filteredTranscriptGtfPath)
        {
            filteredTranscriptGtfPath = Path.Combine(Path.GetDirectoryName(transcriptGtfPath), Path.GetFileNameWithoutExtension(transcriptGtfPath)) + ".filtered" + Path.GetExtension(transcriptGtfPath);
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "echo \"Removing zero-abundance transcripts from " + transcriptGtfPath + "\"",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " ]]; then " +
                    "grep -v 'FPKM \"0.0000000000\"' " + WrapperUtility.ConvertWindowsPath(transcriptGtfPath) + " > " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) +
                "; fi"
            };
        }

        /// <summary>
        /// Merge multiple transcript models (GTF) into a single one (GTF)
        /// </summary>
        /// <param name="spritzDirectory"></param>
        public static List<string> Cuffmerge(string spritzDirectory, int threads, string geneModelGtfOrGff, string genomeFastaPath, List<string> transcriptGtfPaths, string combinedTranscriptGtfOutputPath)
        {
            string gtfListPath = Path.Combine(Path.GetDirectoryName(combinedTranscriptGtfOutputPath), Path.GetFileNameWithoutExtension(combinedTranscriptGtfOutputPath)) + "_gtflist.txt";
            string mergedGtfPath = Path.Combine(combinedTranscriptGtfOutputPath, CufflinksMergedFilename);
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "cufflinks-2.2.1")),
                "readlink -f \"" + String.Join("\" \"", transcriptGtfPaths.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\" > " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(mergedGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(mergedGtfPath) + " ]]; then " +
                    "./cuffmerge" +
                    " -p " + threads.ToString() +
                    " -o " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) +
                    " -g " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGff) +
                    " -s " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) +
                    " --min-isoform-fraction 0 " +
                    WrapperUtility.ConvertWindowsPath(gtfListPath) +
                "; fi"
            };
        }

        /// <summary>
        /// Converts a GFF formatted gene model to GTF
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="geneModelGffPath"></param>
        /// <param name="geneModelGtfPath"></param>
        public static void GffToGtf(string spritzDirectory, string analysisDirectory, string geneModelGffPath, out string geneModelGtfPath)
        {
            if (!Path.GetExtension(geneModelGffPath).StartsWith(".gff"))
            {
                throw new ArgumentException("Input gene model must be gff formatted to convert to gtf.");
            }

            geneModelGtfPath = geneModelGffPath + ".converted.gtf";

            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "GffToGtf.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "echo \"Converting GFF to GTF: " + geneModelGffPath + " -> " + geneModelGtfPath + "\"",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(geneModelGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(geneModelGtfPath) + " ]]; then " +
                    "cufflinks-2.2.1/gffread " + WrapperUtility.ConvertWindowsPath(geneModelGffPath) + " -T -o " + WrapperUtility.ConvertWindowsPath(geneModelGtfPath) +
                "; fi"
            }).WaitForExit();
        }

        #endregion Public Method
    }
}