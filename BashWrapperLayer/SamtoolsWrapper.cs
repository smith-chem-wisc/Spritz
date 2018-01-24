using System;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class SamtoolsWrapper
    {

        public string Samtools { get; set; }

        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installSamtools.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d samtools-1.6 ]; then",
                "  wget --no-check https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2",
                "  tar -jxvf samtools-1.6.tar.bz2",
                "  rm samtools-1.6.tar.bz2",
                "  cd samtools-1.6",
                "  ./configure --prefix=/usr/local/bin",
                "  make",
                "  make install",
                "  chmod +X samtools",
                "  cp samtools /usr/bin/local",
                "fi"
            });
            return scriptPath;
        }

        public static string GenomeFastaIndexCommand(string binDirectory, string genomeFastaPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + ".fai ]; then " + 
                WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "samtools-1.6", "samtools")) + 
                " faidx " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + "; fi";
        }

        public static string IndexBamCommand(string binDirectory, string bamPath)
        {
            return WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "samtools-1.6", "samtools")) + " index " + WrapperUtility.ConvertWindowsPath(bamPath);
        }
    }
}
