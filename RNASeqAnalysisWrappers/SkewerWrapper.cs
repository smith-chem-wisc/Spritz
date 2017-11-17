using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class SkewerWrapper
    {
        public static void trim(string binDirectory, int qualityFilter, string[] readPaths, out string[] readTrimmedPaths, out string log)
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
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(binDirectory),
                WrapperUtility.convert_windows_path(Path.Combine(binDirectory, "skewer-0.2.2", "skewer")) +
                    " -q " + qualityFilter +
                    " -x " + WrapperUtility.convert_windows_path(Path.Combine(binDirectory, "BBMap", "resources", "adapters.fa")) +
                    " " + WrapperUtility.convert_windows_path(readPaths[0]) +
                    (readPaths.Length > 1 ? " " + WrapperUtility.convert_windows_path(readPaths[1]) : ""),
            }).WaitForExit();
            File.Delete(script_path);
        }

        public static void install(string current_directory)
        {
            string script_path = Path.Combine(current_directory, "downloadInstallSkewer.bash");
            if (Directory.Exists(Path.Combine(current_directory, "BBMap")) && Directory.Exists(Path.Combine(current_directory, "skewer-0.2.2"))) return;
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                "git clone https://github.com/BioInfoTools/BBMap.git", // has adapter sequences in the resources file
                "wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz",
                "tar -xvf 0.2.2.tar.gz",
                "rm 0.2.2.tar.gz",
                "cd skewer-0.2.2",
                "make"
            });
            File.Delete(script_path);
        }
    }
}
