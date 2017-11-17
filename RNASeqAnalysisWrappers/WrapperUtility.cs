using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;

namespace RNASeqAnalysisWrappers
{
    public static class WrapperUtility
    {
        private static Regex driveName = new Regex(@"([A-Z]:)");
        private static Regex forwardSlashes = new Regex(@"(\\)");
        public static string ConvertWindowsPath(string path)
        {
            if (path == null) return null;
            if (path == "") return "";
            if (path.StartsWith("/mnt/")) return path;
            return "/mnt/" + Char.ToLowerInvariant(path[0]) + driveName.Replace(forwardSlashes.Replace(path, "/"), "");
        }

        public static Process RunBashCommand(string command, string arguments)
        {
            Process proc = new Process();
            proc.StartInfo.FileName = @"C:\Windows\System32\bash.exe";
            proc.StartInfo.Arguments = "-c \"" + command + " " + arguments +"\"";
            proc.Start();
            return proc;
        }

        public static Process GenerateAndRunScript(string script_path, List<string> commands)
        {
            using (StreamWriter writer = new StreamWriter(script_path))
            {
                writer.Write(AsciiArt() + "\n");
                foreach (string cmd in commands)
                {
                    writer.Write(cmd + "\n");
                }
            }
            return RunBashCommand(@"bash", ConvertWindowsPath(script_path));
        }

        public static bool CheckBashSetup()
        {
            return File.Exists(@"C:\Windows\System32\bash.exe");
        }

        public static void Install(string currentDirectory)
        {
            List<string> commands = new List<string>
            {
                "echo \"Checking for updates and installing any missing dependencies. Please enter your password for this step:\n\"",
                "sudo apt-get update",
                "sudo apt-get upgrade"
            };

            List<string> aptitudeDependencies = new List<string>
            {
                "gcc", "g++", "make", "python", "samtools", "picard-tools", "gawk", "cmake"
            };

            foreach (string dependency in aptitudeDependencies)
            {
                commands.Add(
                    "if commmand -v " + dependency + " > /dev/null 2>&1 ; then\n" +
                    "  echo found\n" +
                    "else\n" +
                    "  sudo apt-get install " + dependency + "\n" +
                    "fi");
            }

            commands.Add(
                "version=$(java -version 2>&1 | awk -F '\"' '/version/ {print $2}')\n" +
                "if [[ \"$version\" > \"1.5\" ]]; then\n" +
                "  echo found\n" +
                "else\n" +
                "  sudo add-apt-repository ppa:webupd8team/java\n" +
                "  sudo apt-get update\n" +
                "  sudo apt-get install oracle-java8-installer\n" +
                "fi");

            string scriptPath = Path.Combine(currentDirectory, "install_dependencies.bash");
            GenerateAndRunScript(scriptPath, commands).WaitForExit();
            //File.Delete(script_path);
        }

        public static void Decompress(string file)
        {
            RunBashCommand("gunzip", file);
        }

        public static string AsciiArt()
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
    }
}
