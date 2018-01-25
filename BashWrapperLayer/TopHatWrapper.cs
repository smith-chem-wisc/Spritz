using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ToolWrapperLayer
{
    public class TopHatWrapper
    {

        #region Public Properties

        public static string BowtieIndexFilePrefix = "index";

        public static string TophatAcceptedHitsFilename = "accepted_hits.bam";

        public static string TophatAlignmentSummaryFilename = "align_summary.txt";

        public static string TophatDeletionsBEDFilename = "deletions.bed";

        public static string TophatJunctionsBEDFilename = "junctions.bed";

        public static string TophatInsertionsBEDFilename = "insertions.bed";

        #endregion Public Properties

        #region Public Methods

        public static bool BowtieIndexExists(string genomeFasta)
        {
            return File.Exists(Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".1.bt2"));
        }

        public static void GenerateBowtieIndex(string binDirectory, string genomeFasta, out string bowtieIndexPrefix)
        {
            bowtieIndexPrefix = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta));
            if (BowtieIndexExists(genomeFasta))
                return;
            string script_name = Path.Combine(binDirectory, "scripts", "bowtieIndices.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "bowtie2-2.3.4/bowtie2-build" +
                    " " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " " + WrapperUtility.ConvertWindowsPath(bowtieIndexPrefix)
            }).WaitForExit();
        }

        public static void Align(string binDirectory, string bowtieIndexPrefix, int threads, string[] fastqPaths, string geneModelGtfOrGffPath, bool strandSpecific, out string outputDirectory)
        {
            string tempDir = Path.Combine(Path.GetDirectoryName(fastqPaths[0]), "tmpDir");
            outputDirectory = Path.Combine(Path.GetDirectoryName(fastqPaths[0]), Path.GetFileNameWithoutExtension(fastqPaths[0]) + "TophatOut");
            Directory.CreateDirectory(tempDir);
            string script_name = Path.Combine(binDirectory, "scripts", "tophatRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "tophat-2.1.1/tophat2" +
                    " --num-threads " + threads.ToString() +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    " --GTF " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " --tmp-dir " + WrapperUtility.ConvertWindowsPath(tempDir) +
                    (strandSpecific ? " --library-type fr-firststrand" : "") +
                    " " + WrapperUtility.ConvertWindowsPath(bowtieIndexPrefix) +
                    " " + String.Join(",", fastqPaths.Select(x => WrapperUtility.ConvertWindowsPath(x)))
            }).WaitForExit();

            if (Directory.Exists(tempDir))
                Directory.Delete(tempDir);
        }

        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installTophat.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d tophat-2.1.1 ]; then",
                "  wget --no-check https://ccb.jhu.edu/software/tophat/downloads/tophat-2.1.1.Linux_x86_64.tar.gz",
                "  wget --no-check -O bowtie.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip/download",
                "  tar -xvf tophat-2.1.1.Linux_x86_64.tar.gz; unzip bowtie.zip",
                "  rm tophat-2.1.1.Linux_x86_64.tar.gz bowtie.zip",
                "  mv tophat-2.1.1.Linux_x86_64 tophat-2.1.1",
                "  mv bowtie2-2.3.4-linux-x86_64 bowtie2-2.3.4",
                "fi"
            });
            return scriptPath;
        }

        #endregion Public Methods

    }
}
