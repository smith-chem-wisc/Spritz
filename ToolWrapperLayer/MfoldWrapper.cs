using System.Collections.Generic;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Mfold calculates the most probable conformations of a DNA or RNA molecule. Can be used to figure out the most probable single-stranded regions.
    /// </summary>
    public class MfoldWrapper :
        IInstallable
    {
        /// <summary>
        /// Writes an installation script for mfold.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallMfold.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d mfold-3.6 ]; then",
                "  wget --no-check http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz",
                "  tar -xvf mfold-3.6.tar.gz",
                "  rm mfold-3.6.tar.gz",
                "  cd mfold-3.6",
                "  ./configure",
                "  make",
                "  sudo make install",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing cufflinks.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }
    }
}