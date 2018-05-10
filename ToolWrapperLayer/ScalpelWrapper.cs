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
        public string IndelVcfPath { get; private set; }

        public string FilteredIndelVcfPath { get; private set; }

        /// <summary>
        /// Write an installation script for scalpel. Requires cmake.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "installScripts", "installScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
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
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        // Need to filter VCF by FILTER = PASS; there are several reasons they don't accept calls that I trust
        // There's an attribute "ZYG" for zygosity, either "het" or "homo" for heterozygous or homozygous
        public List<string> CallIndels(string spritzDirectory, int threads, string genomeFastaP, string bedPath, string bamPath, string outdir)
        {
            var vcftools = new VcfToolsWrapper();
            IndelVcfPath = Path.Combine(outdir, "variants.indel.vcf");
            var commands = new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(IndelVcfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(IndelVcfPath) + " ]]; then " +
                "scalpel-0.5.3/scalpel-discovery --single " +
                    "--bam " + WrapperUtility.ConvertWindowsPath(bamPath) +
                    " --ref " + WrapperUtility.ConvertWindowsPath(genomeFastaP) +
                    " --bed " + WrapperUtility.ConvertWindowsPath(bedPath) +
                    " --numprocs " + threads.ToString() +
                    " --dir " + WrapperUtility.ConvertWindowsPath(outdir) +
                "; fi",
                vcftools.RemoveAllSnvs(spritzDirectory, IndelVcfPath, false, true) // vcf-concat doesn't keep all INFO header lines, so just dump the INFO from each variant
            };
            FilteredIndelVcfPath = vcftools.VcfWithoutSnvsPath;
            return commands;
        }
    }
}