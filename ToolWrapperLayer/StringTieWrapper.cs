using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    public class StringTieWrapper
        : IInstallable
    {
        public List<string> TranscriptGtfPaths { get; private set; }
        public string MergedGtfPath { get; private set; }

        #region Installation Methods

        /// <summary>
        /// Writes a script for installing cufflinks.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installStringTie.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d stringtie-1.3.4d ]; then",
                "  wget --no-check http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz",
                "  tar -xvf stringtie-1.3.4d.Linux_x86_64.tar.gz",
                "  rm stringtie-1.3.4d.Linux_x86_64.tar.gz",
                "  mv stringtie-1.3.4d.Linux_x86_64 stringtie-1.3.4d",
                "  cp stringtie-1.3.4d/stringtie /usr/local/bin",
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
        /// <param name="outputTranscriptGtfPath"></param>
        public static List<string> AssembleTranscripts(string binDirectory, int threads, string bamPath, string geneModelGtfOrGffPath, Genome genome,
            Strandedness strandSpecific, bool inferStrandSpecificity, out string outputTranscriptGtfPath)
        {
            Strandedness strandedness = strandSpecific;
            if (inferStrandSpecificity)
            {
                BAMProperties bamProperties = new BAMProperties(bamPath, geneModelGtfOrGffPath, genome, 0.8);
                strandedness = bamProperties.Strandedness;
            }

            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksSortCheck");
            outputTranscriptGtfPath = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".gtf");
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bamPath) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " + SamtoolsWrapper.SortBam(binDirectory, bamPath) + "; fi",
                "bam=" +  WrapperUtility.ConvertWindowsPath(bamPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then bam=" + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + "; fi",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(outputTranscriptGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(outputTranscriptGtfPath) + " ]]; then",
                "  echo \"Performing stringtie transcript reconstruction on " + bamPath + "\"",
                "  stringtie $bam " +
                    " -p " + threads.ToString() +
                    " -G " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(outputTranscriptGtfPath) +
                    (strandedness == Strandedness.None ? "" : strandedness == Strandedness.Forward ? "--fr" : "--rf"),
                "fi",
            };
        }

        /// <summary>
        /// Merge multiple transcript models (GTF) into a single one (GTF)
        /// </summary>
        /// <param name="binDirectory"></param>
        public static List<string> MergeTranscriptPredictions(string binDirectory, string geneModelGtfOrGffPath, List<string> transcriptGtfPaths, string combinedTranscriptGtfOutputPath)
        {
            string gtfListPath = Path.Combine(Path.GetDirectoryName(combinedTranscriptGtfOutputPath), Path.GetFileNameWithoutExtension(combinedTranscriptGtfOutputPath)) + "_gtflist.txt";
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "cufflinks-2.2.1")),
                "readlink -f \"" + String.Join("\" \"", transcriptGtfPaths.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\" > " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " ]]; then ",
                "  echo \"Performing stringtie transcript merger on GTF list:" + gtfListPath + "\"",
                "  stringtie --merge " +
                    " -G " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) +
                    " " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "fi"
            };
        }

        /// <summary>
        /// Removes transcripts with zero abundance predictions
        /// </summary>
        /// <returns></returns>
        public static List<string> RemoveZeroAbundanceCufflinksPredictionsCommand(string binDirectory, string transcriptGtfPath, out string filteredTranscriptGtfPath)
        {
            filteredTranscriptGtfPath = Path.Combine(Path.GetDirectoryName(transcriptGtfPath), Path.GetFileNameWithoutExtension(transcriptGtfPath)) + ".filtered" + Path.GetExtension(transcriptGtfPath);
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "echo \"Removing zero-abundance transcripts from " + transcriptGtfPath + "\"",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " ]]; then " +
                    "grep -v 'FPKM \"0.000000\"' " + WrapperUtility.ConvertWindowsPath(transcriptGtfPath) + " > " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) +
                "; fi"
            };
        }

        public void TranscriptReconstruction(string binDirectory, string analysisDirectory, int threads, string geneModelGtfOrGff, Genome genome, 
            bool strandSpecific, bool inferStrandSpecificity, List<string> sortedBamFiles)
        {
            string scriptName = Path.Combine(binDirectory, "scripts", "TranscriptReconstruction.bash");
            List<string> reconstructionCommands = new List<string>();
            TranscriptGtfPaths = new List<string>();
            foreach (string sortedBam in sortedBamFiles)
            {
                reconstructionCommands.AddRange(AssembleTranscripts(binDirectory, threads, sortedBam, geneModelGtfOrGff, genome, strandSpecific ? Strandedness.Forward : Strandedness.None, inferStrandSpecificity, out string stringtieGtfTranscriptGtfPath));
                reconstructionCommands.AddRange(RemoveZeroAbundanceCufflinksPredictionsCommand(binDirectory, stringtieGtfTranscriptGtfPath, out string filteredTranscriptModelGtfPath));
                TranscriptGtfPaths.Add(filteredTranscriptModelGtfPath);
            }
            int uniqueSuffix = 1;
            foreach (string f in TranscriptGtfPaths)
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            MergedGtfPath = Path.Combine(analysisDirectory, "MergedStringtieModel" + uniqueSuffix + ".gtf");
            reconstructionCommands.AddRange(MergeTranscriptPredictions(binDirectory, geneModelGtfOrGff, TranscriptGtfPaths, MergedGtfPath));
            WrapperUtility.GenerateAndRunScript(scriptName, reconstructionCommands).WaitForExit();
        }

        #endregion Public Method
    }
}