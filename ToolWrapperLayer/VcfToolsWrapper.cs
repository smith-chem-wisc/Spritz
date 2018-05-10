using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Bedtools is a commonly used toolkit for manipulating BED files.
    /// </summary>
    public class VcfToolsWrapper :
        IInstallable
    {
        public string VcfWithoutIndelsPath { get; private set; }

        public string VcfWithoutSnvsPath { get; private set; }

        public string VcfConcatenatedPath { get; private set; }

        /// <summary>
        /// Writes an install script for bedtools
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "installScripts", "installVcfTools.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "if [ ! -d vcftools-0.1.15 ]; then",
                "  wget --no-check https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz",
                "  tar -xvf vcftools-0.1.15.tar.gz",
                "  rm vcftools-0.1.15.tar.gz",
                "  cd vcftools-0.1.15",
                "  ./configure",
                "  make",
                "  make install",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing bedtools.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        public string RemoveAllIndels(string spritzDirectory, string vcfPath, bool keepInfo, bool applyFilter)
        {
            VcfWithoutIndelsPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath)) + ".NoIndels.vcf";
            return
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(VcfWithoutIndelsPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(VcfWithoutIndelsPath) + " ]; then " +
                    "vcftools " +
                    " --remove-indels --vcf " + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " --recode " +
                    (applyFilter ? " --remove-filtered-all " : "") +
                    (keepInfo ? " --recode-INFO-all " : "") +
                    " --stdout > " + WrapperUtility.ConvertWindowsPath(VcfWithoutIndelsPath) +
                "; fi";
        }

        public string RemoveAllSnvs(string spritzDirectory, string vcfPath, bool keepInfo, bool applyFilter)
        {
            VcfWithoutSnvsPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath)) + ".NoSnvs.vcf";
            return
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(VcfWithoutSnvsPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(VcfWithoutSnvsPath) + " ]; then " +
                    "vcftools " +
                    " --keep-only-indels --vcf " + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " --recode " +
                    (applyFilter ? " --remove-filtered-all " : "") +
                    (keepInfo ? " --recode-INFO-all" : "") +
                    " --stdout > " + WrapperUtility.ConvertWindowsPath(VcfWithoutSnvsPath) +
                "; fi";
        }

        public string Concatenate(string spritzDirectory, IEnumerable<string> vcfInputs, string outPrefix)
        {
            VcfConcatenatedPath = outPrefix + ".concat.vcf";
            return
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(VcfConcatenatedPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(VcfConcatenatedPath) + " ]; then " +
                    "vcf-concat " + string.Join(" ", vcfInputs.Select(v => WrapperUtility.ConvertWindowsPath(v))) + " > " + WrapperUtility.ConvertWindowsPath(VcfConcatenatedPath) +
                "; fi";
        }
    }
}