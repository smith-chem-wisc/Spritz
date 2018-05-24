using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// HISAT2 is a fast and efficient splice RNA-Seq aligner. It has recently (2015) replaced TopHat2 as a low-RAM spliced aligner of choice.
    /// </summary>
    public class HISAT2Wrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes an installation script for HISAT2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallHisat2.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d hisat2-2.1.0 ]; then",
                "  wget --no-check ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip",
                "  unzip hisat2-2.1.0-Linux_x86_64.zip",
                "  rm hisat2-2.1.0-Linux_x86_64.zip",
                "  mv hisat2-2.1.0-Linux_x86_64 hisat2-2.1.0",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing hisat2.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods
        #region Alignment Methods

        public static bool IndexExists(string genomeFasta)
        {
            return File.Exists(Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".1.ht2"));
        }

        public static void GenerateIndex(string spritzDirectory, string analysisDirectory, string genomeFasta, out string IndexPrefix)
        {
            IndexPrefix = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta));
            if (IndexExists(genomeFasta))
            {
                return;
            }
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory,"Hisat2Build.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "hisat2-2.1.0/hisat2-build" + 
                    " " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " " + WrapperUtility.ConvertWindowsPath(IndexPrefix)
            }).WaitForExit();
        }

        public static void Align(string spritzDirectory, string analysisDirectory, string IndexPrefix, string[] fastqPaths, out string outputDirectory)
        {
            outputDirectory = Path.Combine(Path.GetDirectoryName(fastqPaths[0]), Path.GetFileNameWithoutExtension(fastqPaths[0]) + "Hisat2Out.sam");
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "Hisat2Align.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "hisat2-2.1.0/hisat2 -q -x" +
                " " + WrapperUtility.ConvertWindowsPath(IndexPrefix) +
                " -U " + System.String.Join(",", fastqPaths.Select(x => WrapperUtility.ConvertWindowsPath(x))) +
                " -S " + outputDirectory,

            }).WaitForExit();
        }
        #endregion Alignment Methods
    }
}