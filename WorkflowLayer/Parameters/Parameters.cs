using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace WorkflowLayer
{
    public class Parameters
    {
        public Parameters()
        {

        }

        public string Bin { get; set; }
        public string AnalysisDirectory { get; set; }
        public string Reference { get; set; }
        public int Threads { get; set; }
        public List<string[]> Fastq { get; set; }
        public bool StrandSpecific { get; set; }
        public bool InferStrandSpecificity { get; set; }
        public bool OverWriteStarAlignment { get; set; }
        public string GenomeStarIndexDirectory { get; set; }
        public string ReorderedFasta { get; set; }
        public string ProteinFasta { get; set; }
        public string GeneModelGtfOrGff { get; set; }
        public string EnsemblKnownSitesPath { get; set; }
        public List<string> FirstPassSpliceJunctions { get; set; }
        public string SecondPassGenomeDirectory { get; set; }
        public List<string> SortedBamFiles { get; set; }
        public List<string> DedupedBamFiles { get; set; }
        public List<string> ChimericSamFiles { get; set; }
        public List<string> ChimericJunctionFiles { get; set; }
        public bool UseReadSubset { get; set; }
        public int ReadSubset { get; set; }
    }
}
