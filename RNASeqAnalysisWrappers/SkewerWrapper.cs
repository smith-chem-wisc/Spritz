using System.Collections.Generic;
using System.IO;

namespace RNASeqAnalysisWrappers
{
    public class SkewerWrapper
    {
        public void trim(string adapter_sequences_file, int quality_filter, string reads1, string reads2 = "")
        {
            WrapperUtility.run_basic_command(@"skewer", @"-q " + quality_filter + @" -x " + adapter_sequences_file + " " + reads1 + " " + reads2);
        }

        public void download_and_install(string current_directory)
        {
            string script_path = Path.Combine(current_directory, "downloadInstallSkewer.bash");
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
