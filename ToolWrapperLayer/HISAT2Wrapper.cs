using System;
using System.Collections.Generic;
using System.IO;

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
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installHisat2.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
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
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

    }
}
