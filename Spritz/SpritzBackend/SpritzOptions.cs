using System;
using System.IO;
using System.Reflection;

namespace SpritzBackend
{
    public class SpritzOptions
    {
        public string AnalysisDirectory { get; set; }
        public string Fastq1 { get; set; }
        public string Fastq2 { get; set; }
        public string Fastq1SingleEnd { get; set; }
        public string SraAccession { get; set; }
        public string SraAccessionSingleEnd { get; set; }
        public int Threads { get; set; }
        public string Reference { get; set; }
        public bool AnalyzeVariants { get; set; }
        public bool AnalyzeIsoforms { get; set; }
        public bool Quantify { get; set; }
        public bool AvailableReferences { get; set; }
        public bool AnalysisSetup { get; set; }

        public static string DefaultAnalysisDirectory()
        {
            string defaultDirectory = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "results");
            if (!RunnerEngine.IsDirectoryWritable(defaultDirectory))
            {
                defaultDirectory = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "Spritz", "output");
            }
            return defaultDirectory;
        }
    }
}