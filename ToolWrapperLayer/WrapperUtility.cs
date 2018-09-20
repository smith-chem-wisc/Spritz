using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;

namespace ToolWrapperLayer
{
    public static class WrapperUtility
    {
        private static Regex driveName = new Regex(@"([A-Z]:)");
        private static Regex forwardSlashes = new Regex(@"(\\)");

        /// <summary>
        /// Checks if ubuntu bash has been installed
        /// </summary>
        /// <returns></returns>
        public static bool CheckBashSetup()
        {
            string ubuntu = Environment.ExpandEnvironmentVariables(@"%USERPROFILE%\AppData\Local\Microsoft\WindowsApps\ubuntu.exe");
            return File.Exists(ubuntu);
        }

        /// <summary>
        /// Returns true if all tools are set up
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static bool CheckToolSetup(string spritzDirectory)
        {
            return GetToolSetupChecks(spritzDirectory).All(x => x.Item2); 
        }

        /// <summary>
        /// Gets tuples with name of tool and whether the tool is set up
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static List<Tuple<string, bool>> GetToolSetupChecks(string spritzDirectory)
        {
            List<Tuple<string, bool>> setup = new List<Tuple<string, bool>>();
            setup.Add(new Tuple<string, bool>("bedops", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "bedops"))));
            setup.Add(new Tuple<string, bool>("bedtools", File.Exists(Path.Combine(spritzDirectory, "Tools", "bedtools2", "bin", "bedtools"))));
            setup.Add(new Tuple<string, bool>("cufflinks", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "cufflinks-2.2.1"))));
            string gatk = Path.Combine(spritzDirectory, "Tools", "gatk");
            setup.Add(new Tuple<string, bool>("gatk", Directory.Exists(gatk) && Directory.GetFiles(gatk, "gatk*local.jar").Length > 0));
            setup.Add(new Tuple<string, bool>("gatk", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "ChromosomeMappings"))));
            setup.Add(new Tuple<string, bool>("hisat2", File.Exists(Path.Combine(spritzDirectory, "Tools", "hisat2-2.1.0", "hisat2"))));
            setup.Add(new Tuple<string, bool>("mfold", File.Exists(Path.Combine(spritzDirectory, "Tools", "mfold-3.6", "scripts", "mfold"))));
            setup.Add(new Tuple<string, bool>("rsem", File.Exists(Path.Combine(spritzDirectory, "Tools", "RSEM-1.3.0", "rsem-prepare-reference"))));
            setup.Add(new Tuple<string, bool>("rsem", File.Exists(Path.Combine(spritzDirectory, "Tools", "RSEM-1.3.0", "rsem-calculate-expression"))));
            setup.Add(new Tuple<string, bool>("rseqc", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "RSeQC-2.6.4"))));
            setup.Add(new Tuple<string, bool>("samtools", File.Exists(Path.Combine(spritzDirectory, "Tools", "samtools-" + new SamtoolsWrapper().SamtoolsVersion, "samtools"))));
            setup.Add(new Tuple<string, bool>("scalpel", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "scalpel-0.5.3"))));
            setup.Add(new Tuple<string, bool>("scalpel", new ScalpelWrapper().CheckInstallation(spritzDirectory)));
            setup.Add(new Tuple<string, bool>("skewer", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "skewer-0.2.2"))));
            setup.Add(new Tuple<string, bool>("skewer", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "BBMap"))));
            setup.Add(new Tuple<string, bool>("slncky", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "slncky"))));
            setup.Add(new Tuple<string, bool>("slncky", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "slncky", "annotations"))));
            bool toolsFolderExists = Directory.Exists(Path.Combine(spritzDirectory, "Tools"));
            setup.Add(new Tuple<string, bool>("slncky", toolsFolderExists && Directory.GetDirectories(Path.Combine(spritzDirectory, "Tools"), "lastz*").Length > 0));
            setup.Add(new Tuple<string, bool>("snpeff", File.Exists(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.jar"))));
            setup.Add(new Tuple<string, bool>("sratoolkit", toolsFolderExists && Directory.GetDirectories(Path.Combine(spritzDirectory, "Tools"), "sratoolkit*").Length > 0));
            setup.Add(new Tuple<string, bool>("sratoolkit", toolsFolderExists && Directory.GetFiles(
                Directory.GetDirectories(
                    Directory.GetDirectories(Path.Combine(spritzDirectory, "Tools"), "sratoolkit*")[0], "bin")[0], "fastq-dump").Length > 0));
            setup.Add(new Tuple<string, bool>("star", Directory.Exists(Path.Combine(spritzDirectory, "Tools", "STAR-" + STARWrapper.STARVersion))));
            setup.Add(new Tuple<string, bool>("star-fusion", File.Exists(Path.Combine(spritzDirectory, "Tools", "STAR-Fusion-v1.4.0", "STAR-Fusion"))));
            //setup.Add(new Tuple<string, bool>("trinity", Directory.GetDirectories(Path.Combine(spritzDirectory, "Tools"), "trinity*").Length > 0));
            return setup;
        }

        /// <summary>
        /// Converts a Windows-formatted path to UNIX-formatted path
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static string ConvertWindowsPath(string path)
        {
            string pathTrim = path.Trim('"');
            if (pathTrim == null) return null;
            if (pathTrim == "") return "";
            if (pathTrim.StartsWith("/mnt/")) return path;
            return "\"/mnt/" + char.ToLowerInvariant(pathTrim[0]) + driveName.Replace(forwardSlashes.Replace(pathTrim, "/"), "") + "\"";
        }

        /// <summary>
        /// Runs a single bash command from the Ubuntu bash commandline
        /// </summary>
        /// <param name="command"></param>
        /// <param name="arguments"></param>
        /// <returns></returns>
        public static Process RunBashCommand(string command, string arguments)
        {
            Process proc = new Process();
            proc.StartInfo.FileName = Environment.ExpandEnvironmentVariables(@"%USERPROFILE%\AppData\Local\Microsoft\WindowsApps\ubuntu.exe");
            proc.StartInfo.Arguments = $"-c {command} {arguments.Replace(" ", "\\ ")}"; // escape the spaces
            proc.Start();
            return proc;
        }

        /// <summary>
        /// Generates a bash script file
        /// </summary>
        /// <param name="scriptPath"></param>
        /// <param name="commands"></param>
        public static void GenerateScript(string scriptPath, List<string> commands)
        {
            Directory.CreateDirectory(Path.GetDirectoryName(scriptPath));
            using (StreamWriter writer = new StreamWriter(scriptPath))
            {
                writer.Write(SpritzArt() + "\n");
                foreach (string cmd in commands)
                {
                    writer.Write(cmd + "\n");
                }
            }
        }

        /// <summary>
        /// Generates a script file and runs it in Ubuntu bash commandline
        /// </summary>
        /// <param name="scriptPath"></param>
        /// <param name="commands"></param>
        /// <returns></returns>
        public static Process GenerateAndRunScript(string scriptPath, List<string> commands)
        {
            GenerateScript(scriptPath, commands);
            return RunBashCommand(@"bash", ConvertWindowsPath(scriptPath));
        }

        /// <summary>
        /// Gets command to change to the Spritz tools directory
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static string ChangeToToolsDirectoryCommand(string spritzDirectory)
        {
            return "cd " + ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools"));
        }

        /// <summary>
        /// Gets path for an analysis script
        /// </summary>
        /// <param name="analysisDirectory"></param>
        /// <param name="scriptName"></param>
        /// <returns></returns>
        public static string GetAnalysisScriptPath(string analysisDirectory, string scriptName)
        {
            return Path.Combine(analysisDirectory, "Scripts", scriptName);
        }

        /// <summary>
        /// Gets path for an installation script
        /// </summary>
        /// <param name="analysisDirectory"></param>
        /// <param name="scriptName"></param>
        /// <returns></returns>
        public static string GetInstallationScriptPath(string spritzDirectory, string scriptName)
        {
            return Path.Combine(spritzDirectory, "Scripts", "InstallScripts", scriptName);
        }

        /// <summary>
        /// Generates a bash command to check that a file has been closed
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static string EnsureClosedFileCommands(string path)
        {
            return $"exec 3<> {ConvertWindowsPath(path)}; exec 3>&-";
        }

        /// <summary>
        /// Generates some nice art for the original name of this program
        /// </summary>
        /// <returns></returns>
        private static string OldAsciiArt()
        {
            return
                "echo \"" + @"__________                __                _____                    " + "\"\n" +
                "echo \"" + @"\______   \_______  _____/  |_  ____  _____/ ____\___________  _____ " + "\"\n" +
                "echo \"" + @" |     ___/\_  __ \/  _ \   __\/ __ \/  _ \   __\/  _ \_  __ \/     \\" + "\"\n" +
                "echo \"" + @" |    |     |  | \(  <_> )  | \  ___(  <_> )  | (  <_> )  | \/  Y Y  \\" + "\"\n" +
                "echo \"" + @" |____|     |__|   \____/|__|  \___  >____/|__|  \____/|__|  |__|_|  /" + "\"\n" +
                "echo \"" + @"                                   \/                              \/" + "\"\n" +
                "echo \"" + @"________          __        ___.                                     " + "\"\n" +
                "echo \"" + @"\______ \ _____ _/  |______ \_ |__ _____    ______ ____              " + "\"\n" +
                "echo \"" + @" |    |  \\\\__  \\\\   __\__  \ | __ \\\\__  \  /  ___// __ \             " + "\"\n" +
                "echo \"" + @" |    \`   \/ __ \|  |  / __ \| \_\ \/ __ \_\___ \\\\  ___/             " + "\"\n" +
                "echo \"" + @"/_______  (____  /__| (____  /___  (____  /____  >\___  >            " + "\"\n" +
                "echo \"" + @"        \/     \/          \/    \/     \/     \/     \/             " + "\"\n" +
                "echo \"" + @"___________              .__                                         " + "\"\n" +
                "echo \"" + @"\_   _____/ ____    ____ |__| ____   ____                            " + "\"\n" +
                "echo \"" + @" |    __)_ /    \  / ___\|  |/    \_/ __ \                           " + "\"\n" +
                "echo \"" + @" |        \   |  \/ /_/  >  |   |  \  ___/                           " + "\"\n" +
                "echo \"" + @"/_______  /___|  /\___  /|__|___|  /\___  >                          " + "\"\n" +
                "echo \"" + @"        \/     \//_____/         \/     \/                            " + "\"\n";
        }

        /// <summary>
        /// Generates some nice art for Spritz
        /// </summary>
        /// <returns></returns>
        private static string SpritzArt()
        {
            return
                "echo\n" +
                "echo \"" + @"           _                               " + "\"\n" +
                "echo \"" + @"         /' \`\                          /' " + "\"\n" +
                "echo \"" + @"       /'   ._)                     --/'-- " + "\"\n" +
                "echo \"" + @"      (____    ____     ____     O  /'____ " + "\"\n" +
                "echo \"" + @"           ) /'    )--)'    )--/' /' '  _/'" + "\"\n" +
                "echo \"" + @"         /'/'    /' /'       /' /'   _/'   " + "\"\n" +
                "echo \"" + @"(_____,/'/(___,/' /'        (__(___/'__,   " + "\"\n" +
                "echo \"" + @"       /'                                  " + "\"\n" +
                "echo \"" + @"     /'                                    " + "\"\n" +
                "echo \"" + @"   /'                                      " + "\"\n" +
                "echo\n";
        }
    }
}