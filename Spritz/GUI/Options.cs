using System;
using System.IO;
using System.Reflection;
using System.Collections.Generic;

namespace SpritzGUI
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
        public string SnpEff { get; set; }
        public string Organism { get; set; }
        public List<string> Test { get; set; } = new List<string>();

        public Options(int dockerThreads)
        {
            Threads = dockerThreads;
        }
    }
}