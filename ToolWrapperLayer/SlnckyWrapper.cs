using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Slncky is a tool (2016) for annotation of likely lncRNAs based on lack of homology with coding transcripts.
    /// </summary>
    public class SlnckyWrapper :
        IInstallable
    {
        #region Private Fields

        private static string SlnckyAnnotationsLocation = @"https://www.dropbox.com/s/pq7wsjx61sp8ghm/annotations.tar.gz";

        #endregion Private Fields

        #region Installation Methods

        /// <summary>
        /// Writes an installation script for slncky.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installSlncky.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "git clone https://github.com/slncky/slncky.git",
                "cd slncky",
                "if [ ! -d annotations ]; then wget " + SlnckyAnnotationsLocation + "; fi",
                "if [ ! -d annotations ]; then tar -xvf annotations.tar.gz; fi",
                "if [ ! -d annotations ]; then rm annotations.tar.gz; fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing slncky.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Methods

        public void run()
        {
        }

        #endregion Public Methods
    }
}