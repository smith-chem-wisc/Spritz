using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    public class SkewerWrapper
    {
        public static void Trim(string binDirectory, int threads, int qualityFilter, string[] readPaths, out string[] readTrimmedPaths, out string log)
        {
            log = "";
            readTrimmedPaths = new string[readPaths.Length];
            if (readPaths.Length <= 0) return;

            // Only create paired entry if paired input, and ignore inputs after second index
            bool compressed = Path.GetExtension(readPaths[0]) == ".gz";
            string[] uncompressedReadPaths = compressed ? readPaths.Select(x => Path.Combine(Path.GetDirectoryName(x), Path.GetFileNameWithoutExtension(x))).ToArray() : readPaths;
            for (int i = 0; i < readPaths.Length; i++)
            {
                if (i == 0) readTrimmedPaths[0] = Path.Combine(Path.GetDirectoryName(uncompressedReadPaths[0]), Path.GetFileNameWithoutExtension(uncompressedReadPaths[0]) + "-trimmed" + (uncompressedReadPaths.Length > 1 ? "-pair1" : "") + ".fastq");
                if (i == 1) readTrimmedPaths[1] = Path.Combine(Path.GetDirectoryName(uncompressedReadPaths[0]), Path.GetFileNameWithoutExtension(uncompressedReadPaths[0]) + "-trimmed-pair2.fastq");
            }
            log = Path.Combine(Path.GetDirectoryName(uncompressedReadPaths[0]), Path.GetFileNameWithoutExtension(uncompressedReadPaths[0]) + "-trimmed.log");

            bool alreadyTrimmed = File.Exists(readTrimmedPaths[0]) && (readPaths.Length == 1 || File.Exists(readTrimmedPaths[1]));
            if (alreadyTrimmed) return;

            string script_path = Path.Combine(binDirectory, "scripts", "skewered.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "skewer-0.2.2", "skewer")) +
                    " -q " + qualityFilter +
                    " -o " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(uncompressedReadPaths[0]), Path.GetFileNameWithoutExtension(uncompressedReadPaths[0]))) +
                    " -t " + threads.ToString() +
                    " -x " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "BBMap", "resources", "adapters.fa")) +
                    " " + WrapperUtility.ConvertWindowsPath(readPaths[0]) +
                    (readPaths.Length > 1 ? " " + WrapperUtility.ConvertWindowsPath(readPaths[1]) : ""),
            }).WaitForExit();
        }

        public static void Install(string currentDirectory)
        {
            if (Directory.Exists(Path.Combine(currentDirectory, "BBMap")) && Directory.Exists(Path.Combine(currentDirectory, "skewer-0.2.2")))
                return;
            string scriptPath = Path.Combine(currentDirectory, "scripts", "downloadInstallSkewer.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                "git clone https://github.com/BioInfoTools/BBMap.git", // has adapter sequences in the resources file
                "wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz",
                "tar -xvf 0.2.2.tar.gz",
                "rm 0.2.2.tar.gz",
                "cd skewer-0.2.2",
                "make"
            }).WaitForExit();
        }
    }
}
