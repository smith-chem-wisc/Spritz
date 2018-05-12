using System;
using System.Collections.Generic;
using System.Diagnostics;
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
            return File.Exists(@"C:\Windows\System32\bash.exe");
        }

        /// <summary>
        /// Converts a Windows-formatted path to UNIX-formatted path
        /// </summary>
        /// <param name="path"></param>
        /// <returns></returns>
        public static string ConvertWindowsPath(string path)
        {
            if (path == null) return null;
            if (path == "") return "";
            if (path.StartsWith("/mnt/")) return path;
            return "/mnt/" + Char.ToLowerInvariant(path[0]) + driveName.Replace(forwardSlashes.Replace(path, "/"), "");
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
            proc.StartInfo.FileName = @"C:\Windows\System32\bash.exe";
            proc.StartInfo.Arguments = "-c \"" + command + " " + arguments + "\"";
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
            return "exec 3<> " + ConvertWindowsPath(path) + "; exec 3>&-";
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