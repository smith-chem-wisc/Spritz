using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace RNASeqAnalysisWrappers
{
    public class GATKWrapper
    {
        // using BaseRecalibrator
        public static Process base_recalibration()
        {
            return null;
        }

        // using HaplotypeCaller on each BAM file individually
        public static Process variant_calling()
        {
            return null;
        }

        // using GenotypeGVCFs on all VCF files together
        public static Process genotype()
        {
            return null;
        }

        // using VariantRecalibrator, ApplyRecalibration, and then VariantFiltration
        public static Process filter_variants()
        {
            return null;
        }

        public static void install(string current_directory)
        {
            if (Directory.Exists(Path.Combine(current_directory, "GenomeAnalysisTK.jar"))) return;
            string script_path = Path.Combine(current_directory, "install_gatk.sh");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                @"wget https://software.broadinstitute.org/gatk/download/auth?package=GATK",
                "tar -jxvf GenomeAnalysisTK-3.8-0.tar.bz2",
                "rm GenomeAnalysisTK-3.8-0.tar.bz2",
                @"mv GenomeAnalysisTK-*/GenomeAnalysisTK.jar .",
                @"rm -r GenomeAnalysisTK-*",
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
