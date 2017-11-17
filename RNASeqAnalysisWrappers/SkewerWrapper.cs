using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class SkewerWrapper
    {
        public static void Trim(string binDirectory, int threads, int qualityFilter, string[] readPaths, out string[] readTrimmedPaths, out string log)
        {
            log = "";
            readTrimmedPaths = new string[readPaths.Length];
            if (readPaths.Length == 0)
                return;
            readTrimmedPaths[0] = Path.Combine(Path.GetDirectoryName(readPaths[0]), Path.GetFileNameWithoutExtension(readPaths[0]) + "-trimmed" + (readPaths.Length > 1 ? "-pair1" : "") + ".fastq");
            if (readPaths.Length > 1)
                readTrimmedPaths[1] = Path.Combine(Path.GetDirectoryName(readPaths[0]), Path.GetFileNameWithoutExtension(readPaths[0]) + "-trimmed-pair2.fastq");
            log = Path.Combine(Path.GetDirectoryName(readPaths[0]), Path.GetFileNameWithoutExtension(readPaths[0]) + "-trimmed.log");
            bool alreadyTrimmed = File.Exists(readTrimmedPaths[0]) && (readPaths.Length == 1 || File.Exists(readTrimmedPaths[1]));
            if (alreadyTrimmed) return;
            string script_path = Path.Combine(binDirectory, "skewered.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "skewer-0.2.2", "skewer")) +
                    " -q " + qualityFilter +
                    " -t " + threads.ToString() +
                    " -x " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "BBMap", "resources", "adapters.fa")) +
                    " " + WrapperUtility.ConvertWindowsPath(readPaths[0]) +
                    (readPaths.Length > 1 ? " " + WrapperUtility.ConvertWindowsPath(readPaths[1]) : ""),
            }).WaitForExit();
            File.Delete(script_path);
        }

        public static void Install(string currentDirectory)
        {
            string scriptPath = Path.Combine(currentDirectory, "downloadInstallSkewer.bash");
            if (Directory.Exists(Path.Combine(currentDirectory, "BBMap")) && Directory.Exists(Path.Combine(currentDirectory, "skewer-0.2.2"))) return;
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                "git clone https://github.com/BioInfoTools/BBMap.git", // has adapter sequences in the resources file
                "wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz",
                "tar -xvf 0.2.2.tar.gz",
                "rm 0.2.2.tar.gz",
                "cd skewer-0.2.2",
                "make"
            });
            File.Delete(scriptPath);
        }
    }
}
