using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class BEDOPSWrapper
    {
        public static Process gtf2bed(string gtf_path, string current_path)
        {
            string bed_path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".bed");
            if (File.Exists(bed_path)) return null;
            string files = "< " + WrapperUtility.convert_windows_path(gtf_path) + " > " + WrapperUtility.convert_windows_path(bed_path);
            return WrapperUtility.run_basic_command(WrapperUtility.convert_windows_path(Path.Combine(current_path, "bedops", "gtf2bed")), files);
        }

        // see https://github.com/Czh3/NGSTools/blob/master/script/gtf2bed12.sh
        public static void gtf2bed12(string gtf_path, string current_path)
        {
            string tmp_path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".tmp");
            string bed_path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".bed12");
            //if (File.Exists(bed_path)) return;
            string script_path = Path.Combine(current_path, "bed12conversion.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(Path.Combine(current_path, "bedops")),
                "./gtfToGenePred -genePredExt -geneNameAsName2 " + WrapperUtility.convert_windows_path(gtf_path) + " " + WrapperUtility.convert_windows_path(tmp_path),
                "awk '{print $2\"\t\"$4\"\t\"$5\"\t\"$1\"\t0\t\"$3\"\t\"$6\"\t\"$7\"\t0\t\"$8\"\t\"$9\"\t\"$10}' " + WrapperUtility.convert_windows_path(tmp_path) + " > " + WrapperUtility.convert_windows_path(bed_path),
                "rm " + WrapperUtility.convert_windows_path(tmp_path)
            }).WaitForExit();
            File.Delete(script_path);
        }

        public static Process gff2bed(string gff_path, string current_WinPath)
        {
            string bed_path = Path.Combine(Path.GetDirectoryName(gff_path), Path.GetFileNameWithoutExtension(gff_path) + ".bed");
            if (File.Exists(bed_path)) return null;
            string files = "< " + WrapperUtility.convert_windows_path(gff_path) + " > " + WrapperUtility.convert_windows_path(bed_path);
            return WrapperUtility.run_basic_command(WrapperUtility.convert_windows_path(Path.Combine(current_WinPath, "bedops", "gff2bed")), files);
        }

        public static void install(string current_directory)
        {
            if (Directory.Exists(Path.Combine(current_directory, "bedops"))) return;
            string script_path = Path.Combine(current_directory, "install_bedops.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                @"wget https://github.com/bedops/bedops/releases/download/v2.4.29/bedops_linux_x86_64-v2.4.29.tar.bz2",
                "tar -jxvf bedops_linux_x86_64-v2.4.29.tar.bz2",
                "rm bedops_linux_x86_64-v2.4.29.tar.bz2",
                "mv bin bedops",
                "cd bedops",
                "wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred",
                "cd ..",
                "cp bedops/* /usr/local/bin"
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
