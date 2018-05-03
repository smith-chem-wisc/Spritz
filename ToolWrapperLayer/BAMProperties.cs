using Bio.IO.BAM;
using Bio.IO.SAM;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ToolWrapperLayer
{
    public enum RnaSeqProtocol
    {
        SingleEnd,
        PairedEnd,
        Mixture,
    }

    public enum Strandedness
    {
        None,
        Forward,
        Reverse,
    }

    /// <summary>
    /// Gives relevant information about a BAM file
    /// </summary>
    public class BAMProperties
    {
        public BAMProperties(string bamPath, string geneModelPath, Genome genome, double minFractionStrandSpecific)
        {
            CheckProperties(bamPath, geneModelPath, genome, minFractionStrandSpecific);
        }

        public RnaSeqProtocol Protocol { get; private set; }

        public Strandedness Strandedness { get; private set; }

        private double FractionForwardStranded { get; set; } // "forward first strand" library preparation

        private double FractionReverseStranded { get; set; } // "reverse first strand" library preparation

        private double FractionUndetermined { get; set; }

        private Dictionary<string, int> PairedStrandedness { get; } = new Dictionary<string, int>(); // key: readId + mappingStrand + strandFromGene

        private Dictionary<string, int> SingleStrandedness { get; } = new Dictionary<string, int>(); // key: mappingStrand + strandFromGene

        /// <summary>
        /// Given a BAM file, try to guess the RNA-Seq experiment:
		///	1) single-end or pair-end
		///	2) strand_specific or not
		///	3) if it is strand-specific, what's the strand_ness of the protocol
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="bamPath"></param>
        /// <param name="geneModelPath"></param>
        /// <param name="minFractionStrandSpecific"></param>
        /// <returns></returns>
        private void CheckProperties(string bamPath, string geneModelPath, Genome genome, double minFractionStrandSpecific)
        {
            GeneModel gm = new GeneModel(genome, geneModelPath);

            using (var reader = File.OpenRead(bamPath))
            {
                Console.WriteLine("Reading BAM file.");

                // read bam, and filter out reads that are QC failures, unmapped, duplicates, or secondary
                BAMParser bam = new BAMParser();
                var unfilteredReads = bam.Parse(reader).ToList();
                var reads = unfilteredReads.Where(read =>
                    !read.Flag.HasFlag(SAMFlags.QualityCheckFailure) && !read.Flag.HasFlag(SAMFlags.UnmappedQuery)
                        && !read.Flag.HasFlag(SAMFlags.Duplicate) && !read.Flag.HasFlag(SAMFlags.NonPrimeAlignment)).ToList();

                Console.WriteLine("Evaluating reads.");

                Parallel.ForEach(reads, read =>
                {
                    // set the interval contained by this read, and get the gene regions nearby
                    bool isReversed = read.Flag.HasFlag(SAMFlags.QueryOnReverseStrand);
                    Interval readInterval = new Interval(null, read.RName, isReversed ? "-" : "+", read.Pos, read.RefEndPos, null);
                    bool hasNearbyRegion = gm.GenomeForest.Forest.TryGetValue(readInterval.ChromosomeID, out IntervalTree nearbyGeneTree);
                    if (hasNearbyRegion)
                    {
                        List<Interval> nearbyGeneRegions = nearbyGeneTree.Query(readInterval);
                        if (nearbyGeneRegions.Count > 0)
                        {
                            // count up paired-end or single-end read properties
                            string mapStrand = isReversed ? "-" : "+";
                            bool isPaired = read.Flag.HasFlag(SAMFlags.PairedRead);
                            bool isRead1 = read.Flag.HasFlag(SAMFlags.FirstReadInPair);
                            bool isRead2 = read.Flag.HasFlag(SAMFlags.SecondReadInPair);
                            string readId = isRead1 ? "1" : isRead2 ? "2" : null;
                            HashSet<string> strandFromGene = new HashSet<string>(nearbyGeneRegions.Select(x => x.Strand));
                            foreach (string strand in strandFromGene)
                            {
                                Dictionary<string, int> dict = isPaired ? PairedStrandedness : SingleStrandedness;
                                string key = isPaired ?
                                    readId + mapStrand + strand :
                                    mapStrand + strand;
                                lock (dict)
                                {
                                    if (dict.TryGetValue(key, out int count)) { count++; }
                                    else { dict[key] = 1; }
                                }
                            }
                        }
                    }
                });

                // From RSeQC:
                //      Not strand specific:
                // This is PairEnd Data
                // Fraction of reads failed to determine: 0.0172
                // Fraction of reads explained by "1++,1--,2+-,2-+": 0.4903
                // Fraction of reads explained by "1+-,1-+,2++,2--": 0.4925
                //      Strand specific:
                // This is PairEnd Data
                // Fraction of reads failed to determine: 0.0072
                // Fraction of reads explained by "1++,1--,2+-,2-+": 0.9441
                // Fraction of reads explained by "1+-,1-+,2++,2--": 0.0487
                SingleStrandedness.TryGetValue("++", out int sForward1);
                SingleStrandedness.TryGetValue("--", out int sForward2);

                SingleStrandedness.TryGetValue("+-", out int sReverse1);
                SingleStrandedness.TryGetValue("-+", out int sReverse2);

                PairedStrandedness.TryGetValue("1++", out int pForward1);
                PairedStrandedness.TryGetValue("1--", out int pForward2);
                PairedStrandedness.TryGetValue("2+-", out int pForward3);
                PairedStrandedness.TryGetValue("2-+", out int pForward4);

                PairedStrandedness.TryGetValue("1+-", out int pReverse1);
                PairedStrandedness.TryGetValue("1-+", out int pReverse2);
                PairedStrandedness.TryGetValue("2++", out int pReverse3);
                PairedStrandedness.TryGetValue("2--", out int pReverse4);

                if (PairedStrandedness.Count > 0 && SingleStrandedness.Count == 0)
                {
                    Protocol = RnaSeqProtocol.PairedEnd;
                    FractionForwardStranded = (double)(pForward1 + pForward2 + pForward3 + pForward4) / (double)PairedStrandedness.Values.Sum();
                    FractionReverseStranded = (double)(pReverse1 + pReverse2 + pReverse3 + pReverse4) / (double)PairedStrandedness.Values.Sum();
                    FractionUndetermined = 1 - FractionForwardStranded - FractionReverseStranded;
                    if (FractionUndetermined > 0.5)
                    {
                        throw new ArgumentException("A large number of reads failed to determine the standedness of the protocol within " + bamPath);
                    }
                    Strandedness = FractionForwardStranded >= minFractionStrandSpecific ? Strandedness.Forward :
                        FractionReverseStranded >= minFractionStrandSpecific ? Strandedness.Reverse :
                        Strandedness.None;
                }
                else if (SingleStrandedness.Count > 0 && PairedStrandedness.Count == 0)
                {
                    Protocol = RnaSeqProtocol.SingleEnd;
                    FractionForwardStranded = (double)(sForward1 + sForward2) / (double)SingleStrandedness.Values.Sum();
                    FractionReverseStranded = (double)(sReverse1 + sReverse2) / (double)SingleStrandedness.Values.Sum();
                    FractionUndetermined = 1 - FractionForwardStranded - FractionReverseStranded;
                    if (FractionUndetermined > 0.5)
                    {
                        throw new ArgumentException("A large number of reads failed to determine the standedness of the protocol within " + bamPath);
                    }
                    Strandedness = FractionForwardStranded >= minFractionStrandSpecific ? Strandedness.Forward :
                        FractionReverseStranded >= minFractionStrandSpecific ? Strandedness.Reverse :
                        Strandedness.None;
                }
                else
                {
                    Protocol = RnaSeqProtocol.Mixture;
                    Strandedness = Strandedness.None;
                    FractionForwardStranded = (double)(sForward1 + sForward2 + pForward1 + pForward2 + pForward3 + pForward4) / (double)PairedStrandedness.Values.Sum();
                    FractionReverseStranded = (double)(sReverse1 + sReverse2 + pReverse1 + pReverse2 + pReverse3 + pReverse4) / (double)PairedStrandedness.Values.Sum();
                    FractionUndetermined = 1 - FractionForwardStranded - FractionReverseStranded;
                    if (FractionUndetermined > 0.5)
                    {
                        throw new ArgumentException("A large number of reads failed to determine the standedness of the protocol within " + bamPath);
                    }
                    Strandedness = FractionForwardStranded >= minFractionStrandSpecific ? Strandedness.Forward :
                        FractionReverseStranded >= minFractionStrandSpecific ? Strandedness.Reverse :
                        Strandedness.None;
                }
            }
        }

        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            s.AppendLine("Protocol: " + Protocol.ToString());
            s.AppendLine("Strandedness: " + Strandedness.ToString());
            s.AppendLine("Fraction Forward Stranded: " + FractionForwardStranded.ToString());
            s.AppendLine("Fraction Reverse Stranded: " + FractionReverseStranded.ToString());
            s.AppendLine("Fraction Undetermined: " + FractionUndetermined.ToString());
            return s.ToString();
        }
    }
}