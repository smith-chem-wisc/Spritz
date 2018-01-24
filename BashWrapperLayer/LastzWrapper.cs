using System;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class LastzWrapper
    {
        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installLastz.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
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
    }
}
