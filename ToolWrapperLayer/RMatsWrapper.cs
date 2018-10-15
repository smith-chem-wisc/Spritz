using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteogenomics;
using System.IO;

namespace ToolWrapperLayer
{
    public class RMatsWrapper
        : IInstallable
    {
        public string OutputDirectory { get; set; }

        private string RMatsVersion { get; set; } = "4.0.2";

        /// <summary>
        /// Novel junctions between KNOWN exons detected in the RNA-Seq results
        /// Does not contain novel splice sites
        /// </summary>
        public static List<string> GtfNovelEventFilenames { get; } = new List<string>
        {
            "fromGTF.novelEvents.A3SS.txt",
            "fromGTF.novelEvents.A5SS.txt",
            "fromGTF.novelEvents.MXE.txt",
            "fromGTF.novelEvents.RI.txt",
            "fromGTF.novelEvents.SE.txt",
        };

        /// <summary>
        /// Evaluates all possible alternative splicing events derived from GTF
        /// </summary>
        public static List<string> GtfEventFilenames { get; } = new List<string>
        {
            "fromGTF.RI.txt",
            "fromGTF.SE.txt",
            "fromGTF.A3SS.txt",
            "fromGTF.A5SS.txt",
            "fromGTF.MXE.txt",
        };

        /// <summary>
        /// Evaluation of splicing with:
        /// JC uses reads crossing the junction only
        /// JCEC uses reads crossing the junction and reads on the target exon (striped lines in rMATS figure)
        /// </summary>
        public static List<string> RawInputFilenames { get; } = new List<string>
        {
            "JC.raw.input.A3SS.txt",
            "JC.raw.input.A5SS.txt",
            "JC.raw.input.MXE.txt",
            "JC.raw.input.RI.txt",
            "JC.raw.input.SE.txt",
            "JCEC.raw.input.A3SS.txt",
            "JCEC.raw.input.A5SS.txt",
            "JCEC.raw.input.MXE.txt",
            "JCEC.raw.input.RI.txt",
            "JCEC.raw.input.SE.txt",
        };

        /// <summary>
        /// MATS is the statistical analysis, so these have the results of differential splicing analysis
        /// JC uses reads crossing the junction only
        /// JCEC uses reads crossing the junction and reads on the target exon (striped lines in rMATS figure)
        /// </summary>
        public static List<string> MatsResultsFilenames { get; } = new List<string>
        {
            "MXE.MATS.JC.txt",
            "MXE.MATS.JCEC.txt",
            "RI.MATS.JC.txt",
            "RI.MATS.JCEC.txt",
            "SE.MATS.JC.txt",
            "SE.MATS.JCEC.txt",
            "A3SS.MATS.JC.txt",
            "A3SS.MATS.JCEC.txt",
            "A5SS.MATS.JC.txt",
            "A5SS.MATS.JCEC.txt",
        };

        /// <summary>
        /// Dictionary for abbreviations used in filenames. Key: abbreviation. Value: full name.
        /// </summary>
        public static Dictionary<string, string> AlternativeSplicingAbbreviations = new Dictionary<string, string>
        {
            { "SE", "Skipped Exon" },
            { "A5SS", "Alternative 5-Prime Splice Site" },
            { "A3SS", "Alternative 3-Prime Splice Site" },
            { "MXE", "Mutually Exclusive Exons" },
            { "RI", "Retained Intron" },
        };

        /// <summary>
        /// Write an installation script for scalpel. Requires cmake.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string compressedFilename = "rmats." + RMatsVersion + ".tgz";
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d rMATS." + RMatsVersion + " ]; then",
                "  wget --no-check https://sourceforge.net/projects/rnaseq-mats/files/MATS/rMATS." + RMatsVersion + ".tgz/download -O " + compressedFilename,
                "  tar -xzf " + compressedFilename,
                "  rm " + compressedFilename,
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
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallScalpel.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "rm -rf rMATS." + RMatsVersion,
            });
            return scriptPath;
        }

        public string CompareSplicing(string spritzDirectory, string analysisDirectory, int threads, bool isPairedEnd, int readLength,
            string[] bams1, string[] bams2, string customGeneModelGtf)
        {
            int uniqueSuffix = bams1.GetHashCode() ^ bams2.GetHashCode() ^ customGeneModelGtf.GetHashCode();
            OutputDirectory = Path.Combine(analysisDirectory, "rMatsResults" + uniqueSuffix.ToString());

            // write text file with comma-separated list of bam paths for 1 and 2
            // return rMATS command
            // example: python /mnt/e/ProjectsActive/Spritz/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py --nthread 24 
            //   --b1 /mnt/e/ProjectsActive/JurkatProteogenomics/180413/rMatsBam1.txt --b2 /mnt/e/ProjectsActive/JurkatProteogenomics/180413/rMatsBam2.txt 
            //   --gtf /mnt/e/ProjectsActive/JurkatProteogenomics/180413/combined_1-trimmed-pair1Aligned.sortedByCoord.out.filtered.gtf 
            //   --od /mnt/e/ProjectsActive/JurkatProteogenomics/180413/rMatsTest -t paired --readLength 101

            return null;
        }

        public HashSet<string> GetSignificantlyDifferentiallySplicedTranscripts(double maxPValue)
        {
            // read all lines from each of the MATS files, split by tab, GeneID in line[1], PValue in line[18]
            // return all GeneID that had PValue less than max
            return null;
        }
    }
}
