using System;
using System.IO;
using System.Reflection;

namespace SpritzBackend
{
    public class Options : SpritzCmdAppArguments
    {
        public Options()
        {
            AnalysisDirectory = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "results");
            Threads = Environment.ProcessorCount;
            AnalyzeVariants = true;
        }

        public Options(int dockerThreads) : this()
        {
            Threads = dockerThreads;
            if (!RunnerEngine.IsDirectoryWritable(AnalysisDirectory))
            {
                AnalysisDirectory = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "Spritz", "output");
            }
        }
    }
}