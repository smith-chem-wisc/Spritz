using System.Collections.Generic;
using System.IO;
using System.Threading.Tasks;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class InstallFlow
    {

        #region Private Field

        private static List<string> aptitudeDependencies = new List<string>
        {
            // installers
            "gcc",
            "g++",
            "gfortran",
            "make",
            "cmake",
            "build-essential",

            // file compression
            "zlib1g-dev",
            "unzip",

            // bioinformatics
            "samtools",
            "tophat",
            "cufflinks",
            "bedtools",
            "ncbi-blast+",
            "melting",

            // commandline tools
            "gawk",
            "git",
            "python",
            "python-dev",
            "python-setuptools",
            "libpython2.7-dev",
        };

        #endregion Private Field

        #region Public Method

        public static void Run(string binDirectory)
        {
            // get root permissions and update and upgrade the repositories
            List<string> commands = new List<string>
            {
                "echo \"Checking for updates and installing any missing dependencies. Please enter your password for this step:\n\"",
                "sudo apt-get -y update",
                "sudo apt-get -y upgrade",
            };

            // install dependencies from aptitude
            foreach (string dependency in aptitudeDependencies)
            {
                commands.Add(
                    "if commmand -v " + dependency + " > /dev/null 2>&1 ; then\n" +
                    "  echo found\n" +
                    "else\n" +
                    "  sudo apt-get -y install " + dependency + "\n" +
                    "fi");
            }

            // python setup
            commands.Add("sudo easy_install pip");
            commands.Add("sudo pip install --upgrade virtualenv");
            commands.Add("pip install --upgrade pip");
            commands.Add("sudo pip install --upgrade qc bitsets cython bx-python pysam RSeQC numpy"); // for RSeQC

            // java8 setup
            commands.Add(
                "version=$(java -version 2>&1 | awk -F '\"' '/version/ {print $2}')\n" +
                "if [[ \"$version\" > \"1.5\" ]]; then\n" +
                "  echo found\n" +
                "else\n" +
                "  sudo add-apt-repository ppa:webupd8team/java\n" +
                "  sudo apt-get -y update\n" +
                "  sudo apt-get -y install openjdk-8-jdk\n" + // need the JDK for some GATK install, not the JRE found in oracle-java8-installer
                "fi");

            // run some scripts in parallel with root permissions
            string installationLogsDirectory = Path.Combine(binDirectory, "installationLogs");
            Directory.CreateDirectory(installationLogsDirectory);
            List<string> parallelScripts = new List<string>
            {
                // require root permissions
                BEDOPSWrapper.WriteInstallScript(binDirectory),
                LastzWrapper.WriteInstallScript(binDirectory),
                MfoldWrapper.WriteInstallScript(binDirectory),
                GATKWrapper.WriteInstallScript(binDirectory),

                // don't necessarily require root permissions
                RSeQCWrapper.WriteInstallScript(binDirectory),
                ScalpelWrapper.WriteInstallScript(binDirectory),
                SkewerWrapper.WriteInstallScript(binDirectory),
                SlnckyWrapper.WriteInstallScript(binDirectory),
                SnpEffWrapper.WriteInstallScript(binDirectory),
                SRAToolkitWrapper.WriteInstallScript(binDirectory),
                STARWrapper.WriteInstallScript(binDirectory),
                STARFusionWrapper.WriteInstallScript(binDirectory)
            };

            for (int i = 0; i < parallelScripts.Count; i++)
            {
                commands.Add("echo \"Running " + parallelScripts[i] + " in the background.\"");
                commands.Add("bash " + WrapperUtility.ConvertWindowsPath(parallelScripts[i]) + " &> " + 
                    WrapperUtility.ConvertWindowsPath(Path.Combine(installationLogsDirectory, Path.GetFileNameWithoutExtension(parallelScripts[i]) + ".log")) + 
                    " &");
                commands.Add("proc" + i.ToString() + "=$!");
            }

            for (int i = 0; i < parallelScripts.Count; i++)
            {
                commands.Add("wait $proc" + i.ToString());
            }

            // write the and run the installations requiring root permissions
            string scriptPath = Path.Combine(binDirectory, "scripts", "installDependencies.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, commands).WaitForExit();
        }

        #endregion Public Method

    }
}
