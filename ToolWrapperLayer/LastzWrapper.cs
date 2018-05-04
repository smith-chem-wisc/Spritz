using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// LASTZ is a local alignment tool. Used by Slncky in this software.
    /// </summary>
    public class LastzWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes an installation script for lastz.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "installScripts", "installLastz.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "if [ ! -d lastz-1.04.00 ]; then",
                "  wget https://github.com/lastz/lastz/archive/1.04.00.tar.gz",
                "  tar -xvf 1.04.00.tar.gz",
                "  rm 1.04.00.tar.gz",
                "  cd lastz-1.04.00",
                "  make",
                "  chmod +X src/lastz",
                "  sudo cp src/lastz /usr/local/bin",
                "  cd ..",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing lastz.
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