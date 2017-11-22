using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNASeqAnalysisWrappers
{
    public class TopHatWrapper
    {
        public static void generate_index()
        {

        }

        public static int InferInnerDistance

        public static void Align(string binDirectory, bool )
        {
            if (inferStrandSpecificity)
            {
                strandSpecific = RSeQCWrapper.CheckStrandSpecificity(binDirectory, bamPath, geneModelGtfOrGffPath);
            }

            outputDirectory = Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".cufflinksOutput");
            string script_name = Path.Combine(binDirectory, "scripts", "cufflink.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "cufflinks " +
                    " --num-threads " + threads.ToString() +
                    " --GTF-guide " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) +
                    " --output-dir " + WrapperUtility.ConvertWindowsPath(outputDirectory) +
                    (strandSpecific ? "--library-type fr-firststrand" : "") +
                    " " + WrapperUtility.ConvertWindowsPath(bamPath)
            }).WaitForExit();
        }
    }
}
