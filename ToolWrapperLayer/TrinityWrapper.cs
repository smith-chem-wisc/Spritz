using System.Collections.Generic;

namespace ToolWrapperLayer
{
    /// <summary>
    /// LASTZ is a local alignment tool. Used by Slncky in this software.
    /// </summary>
    public class TrinityWrapper :
        IInstallable
    {
        public string TrinityVersion = "2.6.6";

        /// <summary>
        /// Writes an installation script for lastz.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallTrinity.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d trinityrnaseq-Trinity-v" + TrinityVersion + " ]; then",
                "  wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v" + TrinityVersion + ".tar.gz",
                "  tar xvf Trinity-v" + TrinityVersion + ".tar.gz",
                "  rm Trinity-v" + TrinityVersion + ".tar.gz",
                "  cd trinityrnaseq-Trinity-v" + TrinityVersion,
                "  make",
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
    }
}