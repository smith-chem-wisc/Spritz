using System;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class MfoldWrapper
    {
        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installMfold.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
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
    }
}
