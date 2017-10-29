using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
            string script_path = Path.Combine(current_directory, "downloadInstallSkewer.sh");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "git clone https://github.com/BioInfoTools/BBMap.git", // has adapter sequences in the resources file
                "wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz",
                "tar -xvf 0.2.2.tar.gz",
                "rm 0.2.2.tar.gz",
                "cd skewer-0.2.2",
                "make"
            });
        }
    }
}
