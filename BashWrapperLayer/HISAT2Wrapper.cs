using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ToolWrapperLayer
{
    public class HISAT2Wrapper
    {

        public static string WriteInstallScript(string binDirectory)
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
    }
}
