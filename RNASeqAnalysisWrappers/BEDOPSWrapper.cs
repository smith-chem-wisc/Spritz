using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class BEDOPSWrapper
    {
        // Installation takes place in WrapperUtility because it requires sudo access

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

        // see https://gist.github.com/gireeshkbogu/f478ad8495dca56545746cd391615b93
        public static string Gtf2Bed12(string bin, string gtf_path)
        {
            string genePredPath = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".genePred");
            string bed12Path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".bed12");
            string sortedBed12Path = Path.Combine(Path.GetDirectoryName(gtf_path), Path.GetFileNameWithoutExtension(gtf_path) + ".sorted.bed12");
            string scriptPath = Path.Combine(bin, "scripts", "bed12conversion.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "gtfToGenePred " + WrapperUtility.ConvertWindowsPath(gtf_path) + " " + WrapperUtility.ConvertWindowsPath(genePredPath),
                "genePredToBed " + WrapperUtility.ConvertWindowsPath(genePredPath) + " " + WrapperUtility.ConvertWindowsPath(bed12Path),
                "sort -k1,1 -k2,2n " + WrapperUtility.ConvertWindowsPath(bed12Path) + " > " + WrapperUtility.ConvertWindowsPath(sortedBed12Path),
            }).WaitForExit();
            return sortedBed12Path;
        }
    }
}
