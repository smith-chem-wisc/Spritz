using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class NanonetWrapper :
        IInstallable
    {
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallNanonet.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d nanonet ]; then",
                "  git clone https://github.com/nanoporetech/nanonet.git",
                "  cd nanonet",
                "  python setup.py install",
                "fi"
            });
            return scriptPath;
        }

        public string WriteRemoveScript(string spritzDirectory)
        {
            return "";
        }
    }
}