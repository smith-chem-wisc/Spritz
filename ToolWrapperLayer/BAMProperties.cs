using Bio;
using Bio.IO.BAM;
using Bio.IO.Gff;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Bio.IO.SAM;

namespace ToolWrapperLayer
{
    public enum RnaSeqProtocol
    {
        SingleEnd,
        PairedEnd,
        Mixture,
    }

    /// <summary>
    /// Gives relevant information about a BAM file
    /// </summary>
    public class BAMProperties
    {
        public BAMProperties(string bamPath, string geneModelPath, double minFractionStrandSpecific)
        {
            CheckProperties(bamPath, geneModelPath, minFractionStrandSpecific);
        }

        public RnaSeqProtocol Protocol { get; private set; }

        public bool StrandSpecific { get; private set; }

        /// <summary>
        /// Given a BAM file, try to guess the RNA-Seq experiment:
		///	1) single-end or pair-end
		///	2) strand_specific or not
		///	3) if it is strand-specific, what's the strand_ness of the protocol
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="bamPath"></param>
        /// <param name="geneModelPath"></param>
        /// <param name="minFractionStrandSpecific"></param>
        /// <returns></returns>
        private void CheckProperties(string bamPath, string geneModelPath, double minFractionStrandSpecific)
        {
            GeneModel gm = new GeneModel(null, geneModelPath);
            Dictionary<string, int> pairedStrandedness = new Dictionary<string, int>(); // key: readId + mappingStrand + strandFromGene
            Dictionary<string, int> singleStrandedness = new Dictionary<string, int>(); // key: mappingStrand + strandFromGene

            using (var reader = File.OpenRead(bamPath))
            {
                // read bam, and filter out reads that are QC failures, unmapped, duplicates, or secondary
                BAMParser bam = new BAMParser();
                var reads = bam.Parse(reader).Where(read => 
                    !read.Flag.HasFlag(SAMFlags.QualityCheckFailure) && !read.Flag.HasFlag(SAMFlags.UnmappedQuery) 
                        && !read.Flag.HasFlag(SAMFlags.Duplicate) && !read.Flag.HasFlag(SAMFlags.NonPrimeAlignment)).ToList();

                foreach (var read in reads)
                {
                    // set the interval contained by this read, and get the gene regions nearby
                    bool isReversed = read.Flag.HasFlag(SAMFlags.QueryOnReverseStrand);
                    Interval readInterval = new Interval(null, read.QuerySequence.ID, isReversed ? "-" : "+", read.Pos, read.RefEndPos, null);
                    bool hasNearbyRegion = gm.GenomeForest.Forest.TryGetValue(readInterval.ChromosomeID, out IntervalTree nearbyGeneTree);
                    if (!hasNearbyRegion) { continue; }
                    List<Interval> nearbyGeneRegions = nearbyGeneTree.Query(readInterval);
                    if (nearbyGeneRegions.Count == 0) { continue; }
                    HashSet<string> strandFromGene = new HashSet<string>(nearbyGeneRegions.Select(x => x.Strand));

                    // count up paired-end or single-end read properties
                    string mapStrand = isReversed ? "-" : "+";
                    bool isPaired = read.Flag.HasFlag(SAMFlags.PairedRead);
                    bool isRead1 = read.Flag.HasFlag(SAMFlags.FirstReadInPair);
                    bool isRead2 = read.Flag.HasFlag(SAMFlags.SecondReadInPair);
                    string readId = isRead1 ? "1" : isRead2 ? "2" : null;

                    foreach (string strand in strandFromGene)
                    {
                        Dictionary<string, int> dict = isPaired ? pairedStrandedness : singleStrandedness;
                        string key = isPaired ? readId + mapStrand + strand : mapStrand + strand;
                        if (dict.TryGetValue(key, out int count)) { count++; }
                        else { dict[key] = 1; }
                    }
                }

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
                if (pairedStrandedness.Count > 0 && singleStrandedness.Count == 0)
                {
                    Protocol = RnaSeqProtocol.PairedEnd;
                    double spec1 = (double)(pairedStrandedness["1++"] + pairedStrandedness["1--"] + pairedStrandedness["2+-"] + pairedStrandedness["2-+"]) / (double)pairedStrandedness.Values.Sum();
                    double spec2 = (double)(pairedStrandedness["1+-"] + pairedStrandedness["1-+"] + pairedStrandedness["2++"] + pairedStrandedness["2--"]) / (double)pairedStrandedness.Values.Sum();
                    double other = 1 - spec1 - spec2;
                    if (other > 0.5)
                    {
                        throw new ArgumentException("A large number of reads failed to determine the standedness of the protocol within " + bamPath);
                    }
                    StrandSpecific = spec1 < minFractionStrandSpecific || spec2 < minFractionStrandSpecific;
                }
                else if (singleStrandedness.Count > 0 && pairedStrandedness.Count == 0)
                {
                    Protocol = RnaSeqProtocol.SingleEnd;
                    double spec1 = (double)(singleStrandedness["++"] + singleStrandedness["--"]) / (double)singleStrandedness.Values.Sum();
                    double spec2 = (double)(singleStrandedness["+-"] + singleStrandedness["-+"]) / (double)singleStrandedness.Values.Sum();
                    double other = 1 - spec1 - spec2;
                    if (other > 0.5)
                    {
                        throw new ArgumentException("A large number of reads failed to determine the standedness of the protocol within " + bamPath);
                    }
                    StrandSpecific = spec1 < minFractionStrandSpecific || spec2 < minFractionStrandSpecific;
                }
                else
                {
                    Protocol = RnaSeqProtocol.Mixture;
                    StrandSpecific = false;
                }
            }
        }
    }
}