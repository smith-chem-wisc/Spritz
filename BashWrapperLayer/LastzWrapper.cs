using System;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class LastzWrapper
    {
        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installLastz.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d lastz-1.04.00]; then wget https://github.com/lastz/lastz/archive/1.04.00.tar.gz; fi",
                "if [ ! -d lastz-1.04.00]; then tar -xvf 1.04.00.tar.gz; fi",
                "if [ ! -d lastz-1.04.00]; then rm 1.04.00.tar.gz; fi",
                "if [ ! -d lastz-1.04.00]; then cd lastz-1.04.00; fi",
                "if [ ! -d lastz-1.04.00]; then make; fi",
                "if [ ! -d lastz-1.04.00]; then chmod +X src/lastz; fi",
                "if [ ! -d lastz-1.04.00]; then sudo cp src/lastz /usr/local/bin; fi",
                "if [ ! -d lastz-1.04.00]; then cd ..; fi",
            });
            return scriptPath;
        }
    }
}
