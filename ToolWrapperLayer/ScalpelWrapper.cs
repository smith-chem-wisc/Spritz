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

        //public string IndelVcf1IndexedPath { get; private set; }
        public string FilteredIndelVcfPath { get; private set; }

        private string ScalpelLocationCheckFilename { get; set; } = "scalpelLocation.txt";
        private string ScalpelVersion { get; set; } = "0.5.3";

        /// <summary>
        /// Write an installation script for scalpel. Requires cmake.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d scalpel-" + ScalpelVersion + " ]; then",
                "  wget --no-check http://sourceforge.net/projects/scalpel/files/scalpel-" + ScalpelVersion + ".tar.gz; tar zxvf scalpel-" + ScalpelVersion + ".tar.gz; cd scalpel-" + ScalpelVersion + "; make",
                "  cd ..",
                "  rm scalpel-" + ScalpelVersion + ".tar.gz",
                "  readlink -f scalpel-" + ScalpelVersion + " > " + ScalpelLocationCheckFilename,
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
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "RemoveScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "rm -rf scalpel-0.5.3",
            });
            return scriptPath;
        }

        // Need to filter VCF by FILTER = PASS; there are several reasons they don't accept calls that I trust
        // There's an attribute "ZYG" for zygosity, either "het" or "homo" for heterozygous or homozygous
        public List<string> CallIndels(string spritzDirectory, int threads, string genomeFastaP, string bedPath, string bamPath, string outdir)
        {
            CheckInstallation(spritzDirectory);
            var vcftools = new VcfToolsWrapper();
            IndelVcfPath = Path.Combine(outdir, "variants.indel.vcf");
            //IndelVcf1IndexedPath = Path.Combine(outdir, "variants.indel.1index.vcf");
            var commands = new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(IndelVcfPath) + " || ! -s " + WrapperUtility.ConvertWindowsPath(IndelVcfPath) + " ]]; then ",
                "  scalpel-" + ScalpelVersion + "/scalpel-discovery --single " +
                    "--bam " + WrapperUtility.ConvertWindowsPath(bamPath) +
                    " --ref " + WrapperUtility.ConvertWindowsPath(genomeFastaP) +
                    " --bed " + WrapperUtility.ConvertWindowsPath(bedPath) +
                    " --numprocs " + threads.ToString() +
                    " --dir " + WrapperUtility.ConvertWindowsPath(outdir),

                // scalpel uses 0-indexing, where SnpEff uses 1-indexing, so change this output to match snpeff
                //"  awk 'BEGIN{OFS=\"\t\"}{ if (substr($0, 1, 1) != \"#\") $2=++$2; print $0 }' " + WrapperUtility.ConvertWindowsPath(IndelVcfPath) + " > " + WrapperUtility.ConvertWindowsPath(IndelVcf1IndexedPath),
                "fi",

                // vcf-concat doesn't keep all INFO header lines, so just dump the INFO from each variant
                vcftools.RemoveAllSnvs(spritzDirectory, IndelVcfPath, false, true)
            };
            FilteredIndelVcfPath = vcftools.VcfWithoutSnvsPath;
            return commands;
        }

        /// <summary>
        /// Scalpel installation folder cannot be moved and still work. There must be some hard-coded path references.
        /// If needed, delete and reinstall in another location.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns>true if no fix needed, false if fix performed</returns>
        public bool CheckInstallation(string spritzDirectory)
        {
            // Don't go further if installation hasn't been run at all
            string scalpelLocationFile = Path.Combine(spritzDirectory, "Tools", ScalpelLocationCheckFilename);
            if (!File.Exists(scalpelLocationFile))
                return false; 

            // Remove and reinstall if it moved
            string removeScriptPath = WriteRemoveScript(spritzDirectory);
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "CheckScalpelInstallation.bash");
            string expectedLocation = WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "scalpel-" + ScalpelVersion));
            bool isSame = File.ReadAllText(scalpelLocationFile).TrimEnd() == expectedLocation;
            if (!isSame)
            {
                WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
                {
                    WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                    "bash " + WrapperUtility.ConvertWindowsPath(removeScriptPath),
                    "bash " + WrapperUtility.ConvertWindowsPath(WriteInstallScript(spritzDirectory)),
                }).WaitForExit();
            }
            return isSame;
        }
    }
}