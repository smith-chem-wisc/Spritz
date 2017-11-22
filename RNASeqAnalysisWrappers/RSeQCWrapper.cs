using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class RSeQCWrapper
    {
        #region Public Methods

        public static int InferInnerDistance(string binDirectory, string bamPath, string geneModelPath)
        {
            if (Path.GetExtension(geneModelPath) != ".bed")
            {
                geneModelPath = BEDOPSWrapper.GtfOrGff2Bed6(binDirectory, geneModelPath);
            }
            string script_path = Path.Combine(binDirectory, "scripts", "inferInnerDistance.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "python RSeQC-2.6.4/scripts/inner_distance.py" + 
                    " -i " + WrapperUtility.ConvertWindowsPath(bamPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath))) + 
                    " -r " + WrapperUtility.ConvertWindowsPath(geneModelPath)
            }).WaitForExit();

            return 1;
        }

        public static bool CheckStrandSpecificity(string binDirectory, string bamPath, string geneModelPath)
        {
            string outfile = Path.GetFileNameWithoutExtension(bamPath) + ".inferexpt";
            InferExperiment(binDirectory, bamPath, geneModelPath, Path.Combine(binDirectory, outfile));
            string[] lines = File.ReadAllLines(Path.Combine(binDirectory, outfile));
            File.Delete(outfile);
            double fraction_aligned_in_same_direction = double.Parse(lines[lines.Length - 2].Split(':')[1].TrimStart());
            double fraction_aligned_in_other_direction = double.Parse(lines[lines.Length - 1].Split(':')[1].TrimStart());
            return fraction_aligned_in_same_direction / fraction_aligned_in_other_direction < 0.2
                || fraction_aligned_in_same_direction / fraction_aligned_in_other_direction > 0.8;
        }

        public static void Install(string binDirectory)
        {
            if (Directory.Exists(Path.Combine(binDirectory, "RSeQC-2.6.4"))) return;
            string script_path = Path.Combine(binDirectory, "scripts", "install_rseqc.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "wget https://downloads.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz",
                "tar -xvf RSeQC-2.6.4.tar.gz", // infer_experiment.py is in the scripts folder
                "rm RSeQC-2.6.4.tar.gz",
            }).WaitForExit();
            BEDOPSWrapper.Install(binDirectory);
            File.Delete(script_path);
        }

        #endregion Public Methods

        #region Private Methods

        private static void InferExperiment(string binDirectory, string bamFile, string geneModel, string outFile)
        {
            if (Path.GetExtension(geneModel) != ".bed")
            {
                geneModel = BEDOPSWrapper.GtfOrGff2Bed6(binDirectory, geneModel);
            }
            string script_path = Path.Combine(binDirectory, "scripts", "infer_expt.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "python RSeQC-2.6.4/scripts/infer_experiment.py" + 
                    " -r " + WrapperUtility.ConvertWindowsPath(geneModel) + 
                    " -i " + WrapperUtility.ConvertWindowsPath(bamFile) + 
                    " > " + WrapperUtility.ConvertWindowsPath(outFile),
            }).WaitForExit();
        }

        #endregion Private Methods

    }
}
