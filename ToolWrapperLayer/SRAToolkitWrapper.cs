using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// SRAToolkit is used to download sequence read archives from GEO-SRA.
    /// </summary>
    public class SRAToolkitWrapper :
        IInstallable
    {
        private static string AscpDownloadLocation = "http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz";
        private static string DownloadLocation = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz";

        public string[] FastqPaths { get; private set; }
        public string LogPath { get; private set; }

        /// <summary>
        /// Writes an installation script for SRAToolkit.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installSRAToolkit.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if ls sratoolkit* 1> /dev/null 2>&1; then", // if there are files listed matching the pattern sratoolkit*
                "  echo \"Found SRAToolkit.\"",
                "else",
                "  wget " + DownloadLocation,
                "  tar -xvf sratoolkit.current-ubuntu64.tar.gz",
                "  rm sratoolkit.current-ubuntu64.tar.gz",
                "fi"

                //"wget " + AscpDownloadLocation,
                //"tar -xvf aspera-connect*.tar.gz",
                //"rm apera-connect*.tar.gz",
                //"echo \"Installing Aspera ASCP Download Client\"",
                //"sudo sh aspera-connect*.sh",
                //"rm ascp-install-3.5.4.102989-linux-64.sh",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing SRAToolkit.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        public void Fetch(string bin, string sraAccession, string destinationDirectoryPath)
        {
            LogPath = Path.Combine(destinationDirectoryPath, sraAccession + "download.log");
            FastqPaths = Directory.GetFiles(destinationDirectoryPath, sraAccession + "*.fastq");
            if (FastqPaths.Length > 0) // already downloaded
            {
                FastqPaths = FastqPaths.Where(x => x != null && !x.Contains("trimmed") && x.EndsWith(".fastq")).ToArray();
                return;
            };
            string scriptPath = Path.Combine(bin, "scripts", "download" + sraAccession + ".bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "echo \"Downloading " + sraAccession + "\"",
                "cd " + WrapperUtility.ConvertWindowsPath(bin),
                "sratoolkit*/bin/fastq-dump --split-files --outdir \"" + WrapperUtility.ConvertWindowsPath(destinationDirectoryPath) + "\" " +
                    sraAccession + " >> " + WrapperUtility.ConvertWindowsPath(LogPath),
            }).WaitForExit();
            FastqPaths = Directory.GetFiles(destinationDirectoryPath, sraAccession + "*.fastq").ToArray();
        }
    }
}