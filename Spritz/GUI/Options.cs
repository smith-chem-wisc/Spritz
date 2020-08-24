using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;

namespace Spritz
{
    public class Options
    {
        public string AnalysisDirectory { get; set; } = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "output");
        public string Fastq1 { get; set; }
        public string Fastq2 { get; set; }

        //public string ExperimentType { get; set; }
        public string SraAccession { get; set; }

        public int Threads { get; set; } = Environment.ProcessorCount;
        public string Reference { get; set; }

        // new for snakemake
        public string Release { get; set; }

        public string Species { get; set; }
        public string Organism { get; set; }
        public List<string> Test { get; set; } = new List<string>();
        public bool AnalyzeVariants { get; set; } = true;
        public bool AnalyzeIsoforms { get; set; }
        public string SpritzVersion { get; set; }

        public Options(int dockerThreads)
        {
            Threads = dockerThreads;
            try
            {
                string testDirectory = Path.Combine(AnalysisDirectory, $"TestSpritzPermissions{AnalysisDirectory.GetHashCode()}");
                Directory.CreateDirectory(testDirectory);
                Directory.Delete(testDirectory);
            }
            catch (Exception)
            {
                AnalysisDirectory = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "Spritz", "output");
            }
        }
    }
}