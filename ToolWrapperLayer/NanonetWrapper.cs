using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

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
                "  cd nanoneet",
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
