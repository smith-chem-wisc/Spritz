using System.Collections.Generic;
using System.IO;
using System.Linq;
using System;

namespace ToolWrapperLayer
{
    public class SRAToolkitWrapper
    {

        #region Private Fields

        private static string ascpDownloadLocation = "http://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz";
        private static string downloadLocation = "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz";

        #endregion Private Fields

        #region Public Methods

        public static void Fetch(string bin, string sraAccession, string destinationDirectoryPath, out string[] fastqPaths, out string logPath)
        {
            logPath = Path.Combine(destinationDirectoryPath, sraAccession + "download.log");
            fastqPaths = Directory.GetFiles(destinationDirectoryPath, sraAccession + "*.fastq");
            if (fastqPaths.Length > 0) // already downloaded
            {
                fastqPaths = fastqPaths.Where(x => x != null && !x.Contains("trimmed") && x.EndsWith(".fastq")).ToArray();
                return;
            };
            string scriptPath = Path.Combine(bin, "scripts", "download" + sraAccession + ".bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "echo \"Downloading " + sraAccession + "\"",
                "cd " + WrapperUtility.ConvertWindowsPath(bin),
                "sratoolkit*/bin/fastq-dump --split-files --outdir \"" + WrapperUtility.ConvertWindowsPath(destinationDirectoryPath) + "\" " +
                    sraAccession + " > " + WrapperUtility.ConvertWindowsPath(logPath),
            }).WaitForExit();
            fastqPaths = Directory.GetFiles(destinationDirectoryPath, sraAccession + "*.fastq").ToArray();
        }

        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installSRAToolkit.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if ls sratoolkit* 1> /dev/null 2>&1; then", // if there are files listed matching the pattern sratoolkit*
                "  echo \"Found SRAToolkit.\"",
                "else",
                "  wget " + downloadLocation,
                "  tar -xvf sratoolkit.current-ubuntu64.tar.gz",
                "  rm sratoolkit.current-ubuntu64.tar.gz",
                "fi"
                
                //"wget " + ascpDownloadLocation,
                //"tar -xvf aspera-connect*.tar.gz",
                //"rm apera-connect*.tar.gz",
                //"echo \"Installing Aspera ASCP Download Client\"",
                //"sudo sh aspera-connect*.sh",
                //"rm ascp-install-3.5.4.102989-linux-64.sh",
            });
            return scriptPath;
        }

        #endregion Public Methods

    }
}
