using System;
using System.Collections.Generic;
using System.IO;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNASeqAnalysisWrappers
{
    public class BEDOPSWrapper
    {
        public static Process gtf2bed(string gtf_path, string current_path)
        {
            string bed_path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".bed");
            if (File.Exists(bed_path)) return null;
            string files = "< " + WrapperUtility.convert_windows_path(gtf_path) + " > " + WrapperUtility.convert_windows_path(bed_path);
            return WrapperUtility.run_basic_command(WrapperUtility.convert_windows_path(Path.Combine(current_path, "bedops", "gff2bed")), files);
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
            string script_path = Path.Combine(current_directory, "install_bedops.sh");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                @"wget https://github.com/bedops/bedops/releases/download/v2.4.29/bedops_linux_x86_64-v2.4.29.tar.bz2",
                "tar -jxvf bedops_linux_x86_64-v2.4.29.tar.bz2",
                "rm bedops_linux_x86_64-v2.4.29.tar.bz2",
                "mv bin bedops",
                "cp bin/* /usr/local/bin"
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
