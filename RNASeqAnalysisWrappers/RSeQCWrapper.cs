using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace RNASeqAnalysisWrappers
{
    public class RSeQCWrapper
    {
        #region Public Properties

        public static string InnerDistanceRPlotSuffix { get; } = ".inner_distance_plot.r";

        public static string InnerDistanceFrequencyTableSuffix { get; } = ".inner_distance_freq.txt";

        public static string InnerDistanceDistanceTableSuffix { get; } = ".inner_distance.txt";

        #endregion Public Properties

        #region Public Methods

        public static int InferInnerDistance(string binDirectory, string bamPath, string geneModelPath, out string[] outputFiles)
        {
            if (Path.GetExtension(geneModelPath) != ".bed")
            {
                geneModelPath = BEDOPSWrapper.Gtf2Bed12(binDirectory, geneModelPath);
            }
            string script_path = Path.Combine(binDirectory, "scripts", "inferInnerDistance.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "python RSeQC-2.6.4/scripts/inner_distance.py" + 
                    " -i " + WrapperUtility.ConvertWindowsPath(bamPath) + // input
                    " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath))) + // out prefix
                    " -r " + WrapperUtility.ConvertWindowsPath(geneModelPath) // gene model in BED format
            }).WaitForExit();

            outputFiles = new string[]
            {
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceRPlotSuffix,
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceFrequencyTableSuffix,
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceDistanceTableSuffix
            };

            string[] distance_lines = File.ReadAllLines(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceDistanceTableSuffix);
            List<int> distances = new List<int>();
            foreach (string dline in distance_lines)
            {
                if (int.TryParse(dline.Split('\t')[1], out int distance) 
                    && distance < 250 && distance > -250) // default settings for infer_distance
                    distances.Add(distance);
            }
            int averageDistance = (int)Math.Round(distances.Average(), 0);
            return averageDistance;
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
