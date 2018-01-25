using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
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
        /// Writes an install script for bedtools.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installBedtools.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d bedtools2 ]; then",
                "  wget --no-check https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz",
                "  tar -xvf bedtools-2.27.1.tar.gz",
                "  rm bedtools-2.27.1.tar.gz",
                "  cd bedtools2",
                "  make",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing bedtools.
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
