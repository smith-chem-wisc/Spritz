using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ToolWrapperLayer
{
    public class MfoldWrapper
    {
        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installMfold.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d mfold-3.6]; then wget --no-check http://unafold.rna.albany.edu/download/mfold-3.6.tar.gz; fi",
                "if [ ! -d mfold-3.6]; then tar -xvf mfold-3.6.tar.gz; fi",
                "if [ ! -d mfold-3.6]; then rm mfold-3.6.tar.gz; fi",
                "if [ ! -d mfold-3.6]; then cd mfold-3.6; fi",
                "if [ ! -d mfold-3.6]; then ./configure; fi",
                "if [ ! -d mfold-3.6]; then make; fi",
                "if [ ! -d mfold-3.6]; then sudo make install; fi"
            });
            return scriptPath;
        }
    }
}
