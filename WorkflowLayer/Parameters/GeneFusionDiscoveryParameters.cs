using System.Collections.Generic;

namespace WorkflowLayer
{
    public class GeneFusionDiscoveryParameters
        : ISpritzParameters
    {
        public string SpritzDirectory { get; set; }
        public string AnalysisDirectory { get; set; }
        public string Reference { get; set; }
        public int Threads { get; set; } = 1;
        public List<string[]> Fastqs { get; set; }
        public string Organism { get; set; }
        public int MinPeptideLength { get; set; } = 7;
    }
}