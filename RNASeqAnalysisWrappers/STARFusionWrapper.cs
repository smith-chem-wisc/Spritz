using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace RNASeqAnalysisWrappers
{
    public class STARFusionWrapper
    {
        public static void RunStarFusion(string bin_directory, string reference, int threads, string chemericOutJunction, string[] fastq_files, string outdir)
        {
            bool g37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool g38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);

            if (!g37 && !g38) return;

            DownloadFusionPlugNPlay(bin_directory, reference);

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
                    (g37 ? WrapperUtility.ConvertWindowsPath(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) :
                        WrapperUtility.ConvertWindowsPath(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"))) +
                " --chimeric_junction " + WrapperUtility.ConvertWindowsPath(chemericOutJunction);

            string script_name = Path.Combine(bin_directory, "scripts", "star_fusion.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "STAR-Fusion_v1.1.0/STAR-Fusion" + arguments
            }).WaitForExit();
        }

        public static void DownloadFusionPlugNPlay(string binDirectory, string reference)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase) && !Directory.Exists(Path.Combine(binDirectory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play"));
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase) && !Directory.Exists(Path.Combine(binDirectory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"));

            string script_path = Path.Combine(binDirectory, "scripts", "downloadFusionPlugNPlay.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "STAR-Fusion_v1.1.0")),

                downloadGrch37 ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch37 ? "tar -xvf GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch37 ? "rm GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",

                downloadGrch38 ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch38 ? "tar -xvf GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                downloadGrch38 ? "rm GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
            }).WaitForExit();
        }

        public static void Install(string binDirectory)
        {
            bool downloadFusion = !Directory.Exists(Path.Combine(binDirectory, "STAR-Fusion_v1.1.0"));
            string script_path = Path.Combine(binDirectory, "scripts", "install_star.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                downloadFusion ? "wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.1.0/STAR-Fusion_v1.1.0.tar.gz" : "",
                downloadFusion ? "tar -xvf STAR-Fusion_v1.1.0.tar.gz" : "",
                downloadFusion ? "rm STAR-Fusion_v1.1.0.tar.gz" : "",
                downloadFusion ? "cd STAR-Fusion_v1.1.0" : "",
                downloadFusion ? "make" : "",
            }).WaitForExit();
        }
    }
}
