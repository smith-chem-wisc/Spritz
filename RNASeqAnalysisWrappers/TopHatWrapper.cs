using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace RNASeqAnalysisWrappers
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
                "bowtie2-build" +
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
                "tophat" +
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

        #endregion Public Methods

    }
}
