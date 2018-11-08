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
        private static readonly int GapBetweenTranscriptsToMergeTogether = 75;

        public string StringtieVersion { get; set; } = "1.3.4d";
        public List<string> TranscriptGtfPaths { get; private set; } = new List<string>();
        public List<string> FilteredTranscriptGtfPaths { get; private set; } = new List<string>();
        public string MergedGtfPath { get; private set; }
        public string FilteredMergedGtfPath { get; private set; }

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
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(sortedCheckPath) + " ]; then " + SamtoolsWrapper.SortBam(bamPath, threads) + "; fi",
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
            bool strandSpecific, bool inferStrandSpecificity, List<string> sortedBamFiles, bool filterEntriesWithZeroAbundanceStringtieEstimates)
        {
            // transcript reconstruction with stringtie (transcripts and quantities used for lncRNA discovery, etc.)
            List<string> reconstructionCommands = new List<string>();
            foreach (string sortedBam in sortedBamFiles)
            {
                reconstructionCommands.AddRange(AssembleTranscripts(spritzDirectory, threads, sortedBam, geneModelGtfOrGff, genome, strandSpecific ? Strandedness.Forward : Strandedness.None, inferStrandSpecificity, out string stringtieGtfTranscriptGtfPath));
                TranscriptGtfPaths.Add(stringtieGtfTranscriptGtfPath);
            }

            // merge the resultant gene models with the reference (used for sample specific databases)
            int uniqueSuffix = 1;
            foreach (string f in TranscriptGtfPaths)
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            MergedGtfPath = MergedGtfPath = Path.Combine(analysisDirectory, "MergedStringtieModel" + uniqueSuffix + ".gtf");
            reconstructionCommands.AddRange(MergeTranscriptPredictions(spritzDirectory, geneModelGtfOrGff, TranscriptGtfPaths, MergedGtfPath));
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "TranscriptReconstruction.bash"), reconstructionCommands).WaitForExit();

            // filter out the transcripts lacking strand information
            foreach (string gtf in TranscriptGtfPaths)
            {
                string filtered = Path.Combine(Path.GetDirectoryName(gtf), Path.GetFileNameWithoutExtension(gtf) + ".filtered.gtf");
                FilterGtfEntriesWithoutStrand(gtf, filtered, filterEntriesWithZeroAbundanceStringtieEstimates);
                FilteredTranscriptGtfPaths.Add(filtered);
            }
            FilteredMergedGtfPath = Path.Combine(Path.GetDirectoryName(MergedGtfPath), Path.GetFileNameWithoutExtension(MergedGtfPath) + ".filtered.gtf");
            FilterGtfEntriesWithoutStrand(MergedGtfPath, FilteredMergedGtfPath, false); // stringtie merged GTFs no longer have abundance values
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
                "readlink -f \"" + string.Join("\" \"", transcriptGtfPaths.Select(f => WrapperUtility.ConvertWindowsPath(f))) +
                    "\" > " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) + " ]]; then ",
                "  echo \"Performing stringtie transcript merger on GTF list:" + gtfListPath + "\"",
                "  stringtie --merge " +
                    " -G " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(combinedTranscriptGtfOutputPath) +
                    " -g " + GapBetweenTranscriptsToMergeTogether.ToString() +
                    //" -T 0 -F 0" + // filtering is done elsewhere; this really does lead to a lot of bad transcript predictions. Just use the default.
                    //" -f 0.01" + // minimum isoform fraction -- use default in stringtie for now
                    " " + WrapperUtility.ConvertWindowsPath(gtfListPath),
                "fi"
            };
        }

        /// <summary>
        /// Filters GTF or GFF entries that lack strand information
        /// </summary>
        /// <param name="gtfPath"></param>
        /// <param name="gtfOutPath"></param>
        public void FilterGtfEntriesWithoutStrand(string gtfPath, string gtfOutPath, bool filterEntriesWithZeroAbundanceStringtieEstimates)
        {
            var chromFeatures = GeneModel.SimplerParse(gtfPath);
            //if (!File.Exists(gtfOutPath))
            //{
            using (var file = File.Create(gtfOutPath))
            {
                var formatter = new GffFormatter();
                foreach (var chromISeq in chromFeatures)
                {
                    List<MetadataListItem<List<string>>> filteredFeatures = new List<MetadataListItem<List<string>>>();
                    bool isMetadata = chromISeq.Metadata.TryGetValue("features", out object featuresObj);
                    if (isMetadata)
                    {
                        bool okayTranscript = false;
                        var features = featuresObj as List<MetadataListItem<List<string>>>;
                        foreach (var feature in features)
                        {
                            if (!feature.SubItems.TryGetValue("strand", out List<string> strandish)) { continue; }
                            var attributes = GeneModel.SplitAttributes(feature.FreeText);
                            if (feature.Key == "transcript")
                            {
                                bool okayFpkm = !filterEntriesWithZeroAbundanceStringtieEstimates ||
                                    attributes.TryGetValue("FPKM", out string fpkm) && double.TryParse(fpkm, out double fpkmValue) && fpkmValue > 0;
                                bool okayTpm = !filterEntriesWithZeroAbundanceStringtieEstimates ||
                                    attributes.TryGetValue("TPM", out string tpm) && double.TryParse(tpm, out double tpmValue) && tpmValue > 0;
                                okayTranscript = okayFpkm && okayTpm;
                            }
                            if (okayTranscript)
                            {
                                filteredFeatures.Add(feature);
                            }
                        }
                    }
                    chromISeq.Metadata["features"] = filteredFeatures;
                }
                formatter.Format(file, chromFeatures);
            }
            //}
        }
    }
}