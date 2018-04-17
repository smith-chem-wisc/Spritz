using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Samtools is a staple of genomic analysis, allowing indexing, sorting, and inspection of alignment maps.
    /// </summary>
    public class SamtoolsWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes a script to install samtools.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
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
                "  ./configure", // configures install to /usr/local/bin and /usr/local/share
                "  make",
                "  make install",
                "fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing samtools.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Methods

        /// <summary>
        /// Gets a string representing the maximum memory that can be used for samtools sort
        /// </summary>
        /// <returns></returns>
        public static string GetSamtoolsMemoryString()
        {
            var megabytes = Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue());
            megabytes = megabytes > 10000 ? 10000 : megabytes; // this is the max samtools sort can take, apparently
            return megabytes + "M";
        }

        public static string SortBam(string binDirectory, string bamPath)
        {
            return "samtools sort -@ " + Environment.ProcessorCount.ToString() + " -m " + GetSamtoolsMemoryString() +
                " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + " " +
                WrapperUtility.ConvertWindowsPath(bamPath);
        }

        public static string GenomeFastaIndexCommand(string binDirectory, string genomeFastaPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + ".fai ]; then " +
                WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "samtools-1.6", "samtools")) +
                " faidx " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + "; fi";
        }

        public static string IndexBamCommand(string binDirectory, string bamPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(bamPath) + ".bai ]; then " +
                WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "samtools-1.6", "samtools")) + " index " + WrapperUtility.ConvertWindowsPath(bamPath) +
                "; fi";
        }

        #endregion Public Methods
    }
}