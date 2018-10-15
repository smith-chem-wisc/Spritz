using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Samtools is a staple of genomic analysis, allowing indexing, sorting, and inspection of alignment maps.
    /// </summary>
    public class SamtoolsWrapper :
        IInstallable
    {
        public string SamtoolsVersion { get; set; } = "1.8";

        #region Installation Methods

        /// <summary>
        /// Writes a script to install samtools.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallSamtools.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d samtools-" + SamtoolsVersion + " ]; then",
                "  wget --no-check https://github.com/samtools/samtools/releases/download/" + SamtoolsVersion + "/samtools-" + SamtoolsVersion + ".tar.bz2",
                "  tar -jxvf samtools-" + SamtoolsVersion + ".tar.bz2",
                "  rm samtools-" + SamtoolsVersion + ".tar.bz2",
                "  cd samtools-" + SamtoolsVersion + "/htslib-" + SamtoolsVersion,
                "  ./configure", // configures install to /usr/local/bin and /usr/local/share
                "  make",
                "  make install",
                "  cd ..",
                "  ./configure",
                "  make",
                "  make install",
                "fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing samtools.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Methods

        /// <summary>
        /// Gets a string representing 90% of the maximum memory per thread that can be used for samtools sort
        /// </summary>
        /// <returns></returns>
        public static string GetSamtoolsMemoryPerThreadString(int threads)
        {
            int megabytes = (int)(Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue()) * 0.9 / threads);
            megabytes = megabytes > 10000 ? 10000 : megabytes; // this is the max samtools sort can take, apparently
            return megabytes + "M";
        }

        public static string SortBam(string bamPath)
        {
            return "samtools sort -@ " + Environment.ProcessorCount.ToString() + " -m " + GetSamtoolsMemoryPerThreadString(Environment.ProcessorCount) +
                " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(bamPath), Path.GetFileNameWithoutExtension(bamPath) + ".sorted.bam")) + " " +
                WrapperUtility.ConvertWindowsPath(bamPath);
        }

        public static string GenomeFastaIndexCommand(string genomeFastaPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + ".fai ]; then samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + "; fi";
        }

        public static string IndexBamCommand(string bamPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(bamPath) + ".bai ]; then samtools index " + WrapperUtility.ConvertWindowsPath(bamPath) + "; fi";
        }

        public static string GetSequencesFromFasta(string inputFasta, IEnumerable<string> sequenceNames, string outputFasta)
        {
            if (sequenceNames.Any(name => name.Contains(" "))) { throw new ArgumentException("A sequence name query had a space in it; this is not supported" + string.Join(",", sequenceNames) + "."); }
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(outputFasta) + " ]; then samtools faidx " + WrapperUtility.ConvertWindowsPath(inputFasta) +
                " " + string.Join(" ", sequenceNames) + " > " + WrapperUtility.ConvertWindowsPath(outputFasta) + "; fi";
        }

        #endregion Public Methods
    }
}