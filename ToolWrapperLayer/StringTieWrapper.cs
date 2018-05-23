using Bio;
using Bio.IO.Gff;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    public class StringtieWrapper
        : IInstallable
    {
        private static readonly int GapBetweenTranscriptsToMergeTogether = 50;

        public string StringtieVersion { get; set; } = "1.3.4d";
        public List<string> TranscriptGtfPaths { get; private set; } = new List<string>();
        public List<string> FilteredTranscriptGtfPaths { get; private set; } = new List<string>();

        /// <summary>
        /// Writes a script for installing cufflinks.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallStringTie.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d stringtie-" + StringtieVersion + " ]; then",
                "  wget --no-check http://ccb.jhu.edu/software/stringtie/dl/stringtie-" + StringtieVersion + ".Linux_x86_64.tar.gz",
                "  tar -xvf stringtie-" + StringtieVersion + ".Linux_x86_64.tar.gz",
                "  rm stringtie-" + StringtieVersion + ".Linux_x86_64.tar.gz",
                "  mv stringtie-" + StringtieVersion + ".Linux_x86_64 stringtie-" + StringtieVersion,
                "  cp stringtie-" + StringtieVersion + "/stringtie /usr/local/bin",
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
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "rm -rf stringtie-" + StringtieVersion,
            });
            return scriptPath;
        }

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
        /// <param name="outputTranscriptGtfPath"></param>
        public static List<string> AssembleTranscripts(string spritzDirectory, int threads, string bamPath, string geneModelGtfOrGffPath, Genome genome,
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
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bamPath) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " + SamtoolsWrapper.SortBam(spritzDirectory, bamPath) + "; fi",
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
        /// <param name="spritzDirectory"></param>
        public static List<string> MergeTranscriptPredictions(string spritzDirectory, string geneModelGtfOrGffPath, List<string> transcriptGtfPaths,
            string combinedTranscriptGtfOutputPath)
        {
            string gtfListPath = Path.Combine(Path.GetDirectoryName(combinedTranscriptGtfOutputPath), Path.GetFileNameWithoutExtension(combinedTranscriptGtfOutputPath)) +
                "_gtflist.txt";
            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "readlink -f \"" + String.Join("\" \"", transcriptGtfPaths.Select(f => WrapperUtility.ConvertWindowsPath(f))) +
                    "\" > " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " ]]; then ",
                "  echo \"Performing stringtie transcript merger on GTF list:" + gtfListPath + "\"",
                "  stringtie --merge " +
                    " -G " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) +
                    " -g " + GapBetweenTranscriptsToMergeTogether.ToString() +
                    " " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "fi"
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
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "echo \"Removing zero-abundance transcripts from " + transcriptGtfPath + "\"",
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) + " ]]; then " +
                    "grep -v 'FPKM \"0.000000\"' " + WrapperUtility.ConvertWindowsPath(transcriptGtfPath) + " > " + WrapperUtility.ConvertWindowsPath(filteredTranscriptGtfPath) +
                "; fi"
            };
        }

        /// <summary>
        /// Perform transcript reconstruction using stringtie
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="genome"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="sortedBamFiles"></param>
        public void TranscriptReconstruction(string spritzDirectory, string analysisDirectory, int threads, string geneModelGtfOrGff, Genome genome,
            bool strandSpecific, bool inferStrandSpecificity, List<string> sortedBamFiles)
        {
            List<string> reconstructionCommands = new List<string>();
            foreach (string sortedBam in sortedBamFiles)
            {
                reconstructionCommands.AddRange(AssembleTranscripts(spritzDirectory, threads, sortedBam, geneModelGtfOrGff, genome, strandSpecific ? Strandedness.Forward : Strandedness.None, inferStrandSpecificity, out string stringtieGtfTranscriptGtfPath));
                TranscriptGtfPaths.Add(stringtieGtfTranscriptGtfPath);
            }
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "TranscriptReconstruction.bash"), reconstructionCommands).WaitForExit();
            foreach (string gtf in TranscriptGtfPaths)
            {
                string filtered = Path.Combine(Path.GetDirectoryName(gtf), Path.GetFileNameWithoutExtension(gtf) + ".filtered.gtf");
                FilterGtfEntriesWithoutStrand(gtf, filtered);
                FilteredTranscriptGtfPaths.Add(filtered);
            }
        }

        /// <summary>
        /// Filters GTF or GFF entries that lack strand information
        /// </summary>
        /// <param name="gtfPath"></param>
        /// <param name="gtfOutPath"></param>
        public void FilterGtfEntriesWithoutStrand(string gtfPath, string gtfOutPath)
        {
            var chromFeatures = GeneModel.SimplerParse(gtfPath);
            if (!File.Exists(gtfOutPath))
            {
                using (var file = File.Create(gtfOutPath))
                {
                    var formatter = new GffFormatter();
                    foreach (var feature in chromFeatures)
                    {
                        bool isMetadata = feature.Metadata.TryGetValue("features", out object metadata);
                        if (isMetadata)
                        {
                            feature.Metadata["features"] = (metadata as List<MetadataListItem<List<string>>>)
                                .Where(f => f.SubItems.TryGetValue("strand", out List<string> strandish))
                                .ToList();
                            formatter.Format(file, feature);
                        }
                    }
                }
            }
        }
    }
}