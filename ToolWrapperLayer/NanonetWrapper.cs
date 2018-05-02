using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class NanonetWrapper :
        IInstallable
    {
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installNanonet.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d nanonet ]; then",
                "  git clone https://github.com/nanoporetech/nanonet.git",
                "  cd nanonet",
                "  python setup.py install",
                "fi"
            });
            return scriptPath;
        }

        public string WriteRemoveScript(string binDirectory)
        {
            return "";
        }
    }
}