using System.Collections.Generic;
using System.IO;
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
            "liblzma-dev",

            // bioinformatics -- keep this to illustrate that these things are super outdated in aptitude
            //"samtools", // out of date
            //"tophat", // super outdated in aptitude
            //"cufflinks", // super outdated in aptitude
            //"bedtools", // super outdated (2013) in aptitude
            //"ncbi-blast+", // very outdated (2.2 from 2014, instead of 2.7)
            //"melting", // outdated (v4 in aptidude, v5 on main site which was totally rewritten)

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

        public static void Install(string binDirectory)
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
                commands.Add
                (
                    "if commmand -v " + dependency + " > /dev/null 2>&1 ; then\n" +
                    "  echo found\n" +
                    "else\n" +
                    "  sudo apt-get -y install " + dependency + "\n" +
                    "fi"
                );
            }

            // python setup
            commands.Add("sudo easy_install pip");
            commands.Add("sudo pip install --upgrade virtualenv");
            commands.Add("pip install --upgrade pip");
            commands.Add("sudo pip install --upgrade qc bitsets cython bx-python pysam RSeQC numpy"); // for RSeQC

            // java8 setup
            commands.Add
            (
                "version=$(java -version 2>&1 | awk -F '\"' '/version/ {print $2}')\n" +
                "if [[ \"$version\" > \"1.5\" ]]; then\n" +
                "  echo found\n" +
                "else\n" +
                "  sudo add-apt-repository ppa:webupd8team/java\n" +
                "  sudo apt-get -y update\n" +
                "  sudo apt-get -y install openjdk-8-jdk\n" + // need the JDK for some GATK install, not the JRE found in oracle-java8-installer
                "fi"
            );

            // run some scripts in parallel with root permissions
            string installationLogsDirectory = Path.Combine(binDirectory, "scripts", "installLogs");
            Directory.CreateDirectory(installationLogsDirectory);
            List<string> parallelScripts = new List<string>
            {
                // require root permissions
                BEDOPSWrapper.WriteInstallScript(binDirectory),
                BLASTWrapper.WriteInstallScript(binDirectory),
                LastzWrapper.WriteInstallScript(binDirectory),
                MeltingWrapper.WriteInstallScript(binDirectory),
                MfoldWrapper.WriteInstallScript(binDirectory),
                SamtoolsWrapper.WriteInstallScript(binDirectory),

                // don't necessarily require root permissions
                BedtoolsWrapper.WriteInstallScript(binDirectory),
                CufflinksWrapper.WriteInstallScript(binDirectory),
                GATKWrapper.WriteInstallScript(binDirectory),
                HISAT2Wrapper.WriteInstallScript(binDirectory),
                RSeQCWrapper.WriteInstallScript(binDirectory),
                ScalpelWrapper.WriteInstallScript(binDirectory),
                SkewerWrapper.WriteInstallScript(binDirectory),
                SlnckyWrapper.WriteInstallScript(binDirectory),
                SnpEffWrapper.WriteInstallScript(binDirectory),
                SRAToolkitWrapper.WriteInstallScript(binDirectory),
                STARWrapper.WriteInstallScript(binDirectory),
                STARFusionWrapper.WriteInstallScript(binDirectory),
                TopHatWrapper.WriteInstallScript(binDirectory)
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
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installDependencies.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, commands).WaitForExit();
        }

        #endregion Public Method

    }
}
