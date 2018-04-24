using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Scalpel is a tool for accurately calling insertions and deletions.
    /// </summary>
    public class ScalpelWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Write an installation script for scalpel. Requires cmake.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d scalpel-0.5.3 ]; then",
                "  wget --no-check http://sourceforge.net/projects/scalpel/files/scalpel-0.5.3.tar.gz; tar zxvf scalpel-0.5.3.tar.gz; cd scalpel-0.5.3; make",
                "  cd ..",
                "  rm scalpel-0.5.3.tar.gz",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing scalpel.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region Public Method

        // Need to filter VCF by FILTER = PASS; there are several reasons they don't accept calls that I trust
        // There's an attribute "ZYG" for zygosity, either "het" or "homo" for heterozygous or homozygous
        public static List<string> CallIndels(string bin_directory, int threads, string genome_fasta, string bed, string bam, string outdir, out string newVcf)
        {
            newVcf = Path.Combine(outdir, "variants.indel.vcf");
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(newVcf) + " || ! -s " + WrapperUtility.ConvertWindowsPath(newVcf) + " ]]; then " +
                "scalpel-0.5.3/scalpel-discovery --single " +
                    "--bam " + WrapperUtility.ConvertWindowsPath(bam) +
                    " --ref " + WrapperUtility.ConvertWindowsPath(genome_fasta) +
                    " --bed " + WrapperUtility.ConvertWindowsPath(bed) +
                    " --numprocs " + threads.ToString() +
                    " --dir " + WrapperUtility.ConvertWindowsPath(outdir) +
                "; fi",
            };
        }

        #endregion Public Method
    }
}