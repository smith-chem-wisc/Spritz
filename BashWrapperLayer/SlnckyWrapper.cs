using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ToolWrapperLayer
{
    public class SlnckyWrapper
    {
        #region Private Fields

        private static string SlnckyAnnotationsLocation = @"https://www.dropbox.com/s/pq7wsjx61sp8ghm/annotations.tar.gz";

        #endregion Private Fields

        public void run()
        {

        }

        public static void Install(string binDirectory)
        {
            string script_path = Path.Combine(binDirectory, "scripts", "installSlncky.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "git clone https://github.com/slncky/slncky.git",
                "cd slncky",
                "if [ ! -d annotations ]; then wget " + SlnckyAnnotationsLocation + "; fi",
                "if [ ! -d annotations ]; then tar -xvf annotations.tar.gz; fi",
                "if [ ! -d annotations ]; then rm annotations.tar.gz; fi",
            }).WaitForExit();
        }
    }
}
