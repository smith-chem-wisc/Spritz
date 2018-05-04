using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Bedtools is a commonly used toolkit for manipulating BED files.
    /// </summary>
    public class BedtoolsWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes an install script for bedtools 
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "installScripts", "installBedtools.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "if [ ! -d bedtools2 ]; then",
                // (v2.24 is the highest slncky allows)
                "  wget --no-check https://github.com/arq5x/bedtools2/releases/download/v2.24.0/bedtools-2.24.0.tar.gz",
                "  tar -xvf bedtools-2.24.0.tar.gz",
                "  rm bedtools-2.24.0.tar.gz",
                "  cd bedtools2",
                "  make install",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing bedtools.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods
    }
}