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
        private static string SlnckyAnnotationsLocation = @"https://www.dropbox.com/s/pq7wsjx61sp8ghm/annotations.tar.gz";

        public static string CanonicalToLncsSuffix { get; } = ".canonical_to_lncs.txt";

        public static string ClusterInfoSuffix { get; } = ".cluster_info.txt";

        public static string FilteredInfoSuffix { get; } = ".filtered_info.txt";

        public static string LncsBedSuffix { get; } = ".lncs.bed";

        public static string LncsInfoSuffix { get; } = ".lncs.info.txt";

        public static string OrfsSuffix { get; } = ".orfs.txt";

        public static string OrthologsTopSuffix { get; } = ".orthologs.top.txt";

        public static string OrthologsSuffix { get; } = ".orthologs.txt";

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
                "if [ -d annotations ]; then rm annotations.tar.gz; fi",
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

        /// <summary>
        /// Annotate predicted transcripts
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="predictedGeneModelUcscBedPath"></param>
        /// <param name="reference"></param>
        /// <param name="slnckyOutPrefix"></param>
        /// <returns></returns>
        public static List<string> Annotate(string binDirectory, string analysisDirectory, int threads, string predictedGeneModelGtfPath, string reference, string slnckyOutPrefix)
        {
            string sortedBed12Cuffmerge = BEDOPSWrapper.Gtf2Bed12(binDirectory, predictedGeneModelGtfPath);
            string ucscCuffmergeBedPath = EnsemblDownloadsWrapper.ConvertFirstColumnEnsembl2UCSC(binDirectory, reference, sortedBed12Cuffmerge);
            Directory.CreateDirectory(Path.GetDirectoryName(slnckyOutPrefix));
            string ucscReference = reference.Contains("38") ? "hg38" : "hg19";
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "slncky")),
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(slnckyOutPrefix + LncsBedSuffix) + " || ! -s " + WrapperUtility.ConvertWindowsPath(slnckyOutPrefix + LncsBedSuffix) + " ]]; then " +
                    "./slncky.v1.0" +
                    " --threads " + threads.ToString() +
                    " " + WrapperUtility.ConvertWindowsPath(ucscCuffmergeBedPath) +
                    " " + ucscReference +
                    " " + WrapperUtility.ConvertWindowsPath(slnckyOutPrefix) +
                "; fi"
            };
        }

        #endregion Public Methods
    }
}