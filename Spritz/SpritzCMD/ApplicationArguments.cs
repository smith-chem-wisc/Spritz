﻿using System.Collections.Generic;

namespace SpritzCMD
{
    public class ApplicationArguments
    {
        public string AnalysisDirectory { get; set; }
        public string Fastq1 { get; set; }
        public string Fastq2 { get; set; }
        public string Fastq1SingleEnd { get; set; }
        public string SraAccession { get; set; }
        public string SraAccessionSingleEnd { get; set; }
        public int Threads { get; set; }
        public string Reference { get; set; }
        public string Release { get; set; }
        public string Species { get; set; }
        public string Organism { get; set; }
        public bool AnalyzeVariants { get; set; } = true;
        public bool AnalyzeIsoforms { get; set; }
        public bool AvailableReferences { get; set; }
        public bool AnalysisSetup { get; set; }
    }
}