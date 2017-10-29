using System;
using System.Collections.Generic;
using System.IO;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNASeqAnalysisWrappers
{
    public class RSeQCWrapper
    {
        public static Process infer_experiment(string bin_directory, string bam_file, string bed_genes, string out_file)
        {
            return WrapperUtility.run_basic_command("python", WrapperUtility.convert_windows_path(Path.Combine(bin_directory, "RSeQC-2.6.4", "scripts", "infer_experiment.py")) + " -s 100000 -r " + WrapperUtility.convert_windows_path(bed_genes) + " -i " + WrapperUtility.convert_windows_path(bam_file) + " > " + WrapperUtility.convert_windows_path(out_file));
        }

        public static bool check_strand_specificity(string bin_directory, string bam_file, string bed_genes, string current_directory)
        {
            string outfile = Path.GetFileNameWithoutExtension(bam_file) + ".inferexpt";
            infer_experiment(bin_directory, bam_file, bed_genes, Path.Combine(current_directory, outfile))
                .WaitForExit();
            string[] lines = File.ReadAllLines(Path.Combine(current_directory, outfile));
            File.Delete(outfile);
            double fraction_aligned_in_other_direction = double.Parse(lines[lines.Length - 1].Split(':')[1].TrimStart());
            return fraction_aligned_in_other_direction < 0.2;
        }

        public static void install(string current_directory)
        {
            if (Directory.Exists(Path.Combine(current_directory, "RSeQC-2.6.4"))) return;
            string script_path = Path.Combine(current_directory, "install_rseqc.sh");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                "wget https://downloads.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz",
                "tar -xvf RSeQC-2.6.4.tar.gz", // infer_experiment.py is in the scripts folder
                "rm RSeQC-2.6.4.tar.gz", 
                "pip install bx-python pysam"
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
