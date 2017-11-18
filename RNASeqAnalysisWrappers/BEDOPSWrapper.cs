using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class BEDOPSWrapper
    {
        // see https://www.biostars.org/p/206342/ for awk fix
        public static string GtfOrGff2Bed6(string bin, string gtfOrGffPath)
        {
            string extension = Path.GetExtension(gtfOrGffPath);
            string bedPath = Path.Combine(Path.GetDirectoryName(gtfOrGffPath), Path.GetFileNameWithoutExtension(gtfOrGffPath) + ".bed");
            if (!File.Exists(bedPath))
            {
                string scriptPath = Path.Combine(bin, "scripts", "bed6conversion.bash");
                WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
                {
                    "cd " + WrapperUtility.ConvertWindowsPath(bin),
                     (extension == ".gtf" ? "awk '{ if ($0 ~ \"transcript_id\") print $0; else print $0\" transcript_id \\\"\\\";\"; }' " : "cat ") + WrapperUtility.ConvertWindowsPath(gtfOrGffPath) 
                        + " | " + WrapperUtility.ConvertWindowsPath(Path.Combine(bin, "bedops", extension == ".gtf" ? "gtf2bed" : "gff2bed")) + 
                        " - > " + WrapperUtility.ConvertWindowsPath(bedPath),
                }).WaitForExit();
            }
            return bedPath;
        }

        // see https://github.com/Czh3/NGSTools/blob/master/script/gtf2bed12.sh
        public static string Gtf2Bed12(string bin, string gtf_path)
        {
            string tmpPath = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".tmp");
            string bedPath = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".bed12");
            string scriptPath = Path.Combine(bin, "scripts", "bed12conversion.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(bin, "bedops")),
                "./gtfToGenePred -genePredExt -geneNameAsName2 " + WrapperUtility.ConvertWindowsPath(gtf_path) + " " + WrapperUtility.ConvertWindowsPath(tmpPath),
                "awk '{print $2\"\t\"$4\"\t\"$5\"\t\"$1\"\t0\t\"$3\"\t\"$6\"\t\"$7\"\t0\t\"$8\"\t\"$9\"\t\"$10}' " + WrapperUtility.ConvertWindowsPath(tmpPath) + " > " + WrapperUtility.ConvertWindowsPath(bedPath),
                "rm " + WrapperUtility.ConvertWindowsPath(tmpPath)
            }).WaitForExit();
            return bedPath;
        }

        public static void Install(string currentDirectory)
        {
            if (Directory.Exists(Path.Combine(currentDirectory, "bedops"))) return;
            string scriptPath = Path.Combine(currentDirectory, "scripts", "install_bedops.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                @"wget https://github.com/bedops/bedops/releases/download/v2.4.29/bedops_linux_x86_64-v2.4.29.tar.bz2",
                "tar -jxvf bedops_linux_x86_64-v2.4.29.tar.bz2",
                "rm bedops_linux_x86_64-v2.4.29.tar.bz2",
                "mv bin bedops",
                "cd bedops",
                "wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred",
                "cd ..",
                "cp bedops/* /usr/local/bin"
            }).WaitForExit();
        }
    }
}
