using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class RSeQCWrapper
    {
        #region Public Methods

        public static bool CheckStrandSpecificity(string binDirectory, string bamFile, string geneModel)
        {
            string outfile = Path.GetFileNameWithoutExtension(bamFile) + ".inferexpt";
            InferExperiment(binDirectory, bamFile, geneModel, Path.Combine(binDirectory, outfile));
            string[] lines = File.ReadAllLines(Path.Combine(binDirectory, outfile));
            File.Delete(outfile);
            double fraction_aligned_in_same_direction = double.Parse(lines[lines.Length - 2].Split(':')[1].TrimStart());
            double fraction_aligned_in_other_direction = double.Parse(lines[lines.Length - 1].Split(':')[1].TrimStart());
            return fraction_aligned_in_same_direction / fraction_aligned_in_other_direction < 0.2
                || fraction_aligned_in_same_direction / fraction_aligned_in_other_direction > 0.8;
        }

        public static void Install(string currentDirectory)
        {
            if (Directory.Exists(Path.Combine(currentDirectory, "RSeQC-2.6.4"))) return;
            string script_path = Path.Combine(currentDirectory, "install_rseqc.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                "wget https://downloads.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz",
                "tar -xvf RSeQC-2.6.4.tar.gz", // infer_experiment.py is in the scripts folder
                "rm RSeQC-2.6.4.tar.gz",
                "pip install bx-python pysam"
            }).WaitForExit();
            BEDOPSWrapper.Install(currentDirectory);
            File.Delete(script_path);
        }

        #endregion Public Methods

        #region Private Methods

        private static void InferExperiment(string binDirectory, string bamFile, string geneModel, string outFile)
        {
            if (Path.GetExtension(geneModel) != ".bed")
                geneModel = BEDOPSWrapper.GtfOrGff2Bed6(binDirectory, geneModel);
            string script_path = Path.Combine(binDirectory, "infer_expt.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "python RSeQC-2.6.4/scripts/infer_experiment.py -r " + WrapperUtility.ConvertWindowsPath(geneModel) + " -i " + WrapperUtility.ConvertWindowsPath(bamFile) + " > " + WrapperUtility.ConvertWindowsPath(outFile),
            }).WaitForExit();
            File.Delete(script_path);
        }

        #endregion Private Methods

    }
}
