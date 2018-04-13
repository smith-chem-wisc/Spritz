using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// RSeQC is a tookit of quality control scripts.
    /// </summary>
    public class RSeQCWrapper :
        IInstallable
    {
        #region Public Properties

        public static string InnerDistanceRPlotSuffix { get; } = ".inner_distance_plot.r";

        public static string InnerDistanceFrequencyTableSuffix { get; } = ".inner_distance_freq.txt";

        public static string InnerDistanceDistanceTableSuffix { get; } = ".inner_distance.txt";

        #endregion Public Properties

        #region Installation Methods

        /// <summary>
        /// Writes an installation script for RSeQC.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installRSeQC.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d RSeQC-2.6.4 ]; then",
                "  wget https://downloads.sourceforge.net/project/rseqc/RSeQC-2.6.4.tar.gz",
                "  tar -xvf RSeQC-2.6.4.tar.gz", // infer_experiment.py is in the scripts folder
                "  rm RSeQC-2.6.4.tar.gz",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script to remove RSeQC.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="bamPath"></param>
        /// <param name="geneModelPath"></param>
        /// <param name="outputFiles"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        public static int InferInnerDistance(string binDirectory, string bamPath, string geneModelPath, out string[] outputFiles)
        {
            if (Path.GetExtension(geneModelPath) != ".bed")
            {
                geneModelPath = BEDOPSWrapper.Gtf2Bed12(binDirectory, geneModelPath);
            }

            outputFiles = new string[]
            {
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceRPlotSuffix,
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceFrequencyTableSuffix,
                Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath)) + InnerDistanceDistanceTableSuffix
            };

            string script_path = Path.Combine(binDirectory, "scripts", "inferInnerDistance.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "python RSeQC-2.6.4/scripts/inner_distance.py" +
                    " -i " + WrapperUtility.ConvertWindowsPath(bamPath) + // input
                    " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath))) + // out prefix
                    " -r " + WrapperUtility.ConvertWindowsPath(geneModelPath), // gene model in BED format
                WrapperUtility.EnsureClosedFileCommands(outputFiles[0]),
                WrapperUtility.EnsureClosedFileCommands(outputFiles[1]),
                WrapperUtility.EnsureClosedFileCommands(outputFiles[2]),
            }).WaitForExit();

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
    }
}