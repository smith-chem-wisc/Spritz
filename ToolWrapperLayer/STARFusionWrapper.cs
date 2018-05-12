using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// STARFusion is a program for analyzing gene fusion events detected using the STAR aligner.
    /// </summary>
    public class STARFusionWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes a script for installing STAR-Fusion.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallStarFusion.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d STAR-Fusion_v1.1.0 ]; then",
                "  wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.1.0/STAR-Fusion_v1.1.0.tar.gz",
                "  tar -xvf STAR-Fusion_v1.1.0.tar.gz",
                "  rm STAR-Fusion_v1.1.0.tar.gz",
                "  cd STAR-Fusion_v1.1.0",
                "  make",
                "fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing STAR-Fusion.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Methods

        public static void RunStarFusion(string spritzDirectory, string analysisDirectory, string reference, int threads, string chemericOutJunction, string[] fastq_files, string outdir)
        {
            bool g37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool g38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);

            if (!g37 && !g38) return;

            DownloadFusionPlugNPlay(spritzDirectory, analysisDirectory, reference);

            string read_command = fastq_files.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastq_files.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --annotate --examine_coding_effects" +
                (fastq_files.Length > 0 ? " --left_fq " + fastq_files[0] : "") +
                (fastq_files.Length > 1 ? " --right_fq " + fastq_files[1] : "") +
                " --CPU " + threads.ToString() +
                " --output_dir " + WrapperUtility.ConvertWindowsPath(outdir) +
                " --genome_lib_dir " +
                    (g37 ? WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) :
                        WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"))) +
                " --chimeric_junction " + WrapperUtility.ConvertWindowsPath(chemericOutJunction);

            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "StarFusion.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "STAR-Fusion_v1.1.0/STAR-Fusion" + arguments
            }).WaitForExit();
        }

        public static void DownloadFusionPlugNPlay(string spritzDirectory, string analysisDirectory, string reference)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase) && !Directory.Exists(Path.Combine(spritzDirectory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play"));
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase) && !Directory.Exists(Path.Combine(spritzDirectory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"));

            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "DownloadFusionPlugNPlay.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "STAR-Fusion_v1.1.0")),

                downloadGrch37 ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch37 ? "tar -xvf GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch37 ? "rm GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",

                downloadGrch38 ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch38 ? "tar -xvf GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch38 ? "rm GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
            }).WaitForExit();
        }

        #endregion Public Methods
    }
}