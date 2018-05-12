using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// TopHat2 is a spliced aligner for RNA-Seq reads. The authors of the tool recommend using HISAT2 instead; it's more sensitive and efficient.
    /// </summary>
    public class TopHatWrapper :
        IInstallable
    {
        #region Public Properties

        /// <summary>
        /// File prefix for bowtie indices.
        /// </summary>
        public static string BowtieIndexFilePrefix = "index";

        /// <summary>
        /// Output filename for TopHat-aligned reads.
        /// </summary>
        public static string TophatAcceptedHitsFilename = "accepted_hits.bam";

        /// <summary>
        /// Output filename for TopHat summary.
        /// </summary>
        public static string TophatAlignmentSummaryFilename = "align_summary.txt";

        /// <summary>
        /// Output filename for detected deletions.
        /// </summary>
        public static string TophatDeletionsBEDFilename = "deletions.bed";

        /// <summary>
        /// Output filename for detected splice junctions.
        /// </summary>
        public static string TophatJunctionsBEDFilename = "junctions.bed";

        /// <summary>
        /// Output filename for detected insertions.
        /// </summary>
        public static string TophatInsertionsBEDFilename = "insertions.bed";

        #endregion Public Properties

        #region Installation Methods

        /// <summary>
        /// Writes an installation script for TopHat and bowtie2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallTophat.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d tophat-2.1.1 ]; then",
                "  wget --no-check https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz",
                "  wget --no-check -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip/download",
                "  tar -xvf tophat-2.1.1.Linux_x86_64.tar.gz; unzip bowtie.zip",
                "  rm tophat-2.1.1.Linux_x86_64.tar.gz bowtie.zip",
                "  mv tophat-2.1.1.Linux_x86_64 tophat-2.1.1",
                "  mv bowtie2-2.3.4-linux-x86_64 bowtie2-2.3.4",
                "  cp bowtie2-2.3.4/* tophat-2.1.1",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing TopHat and bowtie2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Methods

        /// <summary>
        /// Checks whether the bowtie indices exist.
        /// </summary>
        /// <param name="genomeFasta"></param>
        /// <returns></returns>
        public static bool BowtieIndexExists(string genomeFasta)
        {
            return File.Exists(Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".1.bt2"));
        }

        /// <summary>
        /// Generates indices using bowtie2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="bowtieIndexPrefix"></param>
        public static void GenerateBowtieIndex(string spritzDirectory, string analysisDirectory, string genomeFasta, out string bowtieIndexPrefix)
        {
            bowtieIndexPrefix = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta));
            if (BowtieIndexExists(genomeFasta))
            {
                return;
            }
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "BowtieIndices.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "bowtie2-2.3.4/bowtie2-build" +
                    " " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " " + WrapperUtility.ConvertWindowsPath(bowtieIndexPrefix)
            }).WaitForExit();
        }

        /// <summary>
        /// Aligns reads in fastq files using TopHat2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="bowtieIndexPrefix"></param>
        /// <param name="threads"></param>
        /// <param name="fastqPaths"></param>
        /// <param name="geneModelGtfOrGffPath"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="outputDirectory"></param>
        public static void Align(string spritzDirectory, string analysisDirectory, string bowtieIndexPrefix, int threads, string[] fastqPaths, bool strandSpecific, out string outputDirectory)
        {
            string tempDir = Path.Combine(Path.GetDirectoryName(fastqPaths[0]), "tmpDir");
            outputDirectory = Path.Combine(Path.GetDirectoryName(fastqPaths[0]), Path.GetFileNameWithoutExtension(fastqPaths[0]) + "TophatOut");
            Directory.CreateDirectory(tempDir);
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "TophatRun.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "tophat-2.1.1/tophat2" +
                    " --num-threads " + threads.ToString() +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    //" --GTF " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) + /// this triggers tophat to try building an index
                    " --tmp-dir " + WrapperUtility.ConvertWindowsPath(tempDir) +
                    (strandSpecific ? " --library-type fr-firststrand" : "") +
                    " " + WrapperUtility.ConvertWindowsPath(bowtieIndexPrefix) +
                    " " + String.Join(",", fastqPaths.Select(x => WrapperUtility.ConvertWindowsPath(x))),
                "if [ -d " + WrapperUtility.ConvertWindowsPath(tempDir) + " ]; then rm -r " + WrapperUtility.ConvertWindowsPath(tempDir) + "; fi",
            }).WaitForExit();
        }

        /// <summary>
        /// Gets the Windows-formatted path to the directory containing bowtie2
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static string GetBowtie2DirectoryPath(string spritzDirectory)
        {
            return Path.Combine(spritzDirectory, "tophat-2.1.1");
        }

        #endregion Public Methods
    }
}