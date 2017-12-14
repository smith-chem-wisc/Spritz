using Bio;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class Transcript
    {

        #region Private Properties

        private MetadataListItem<List<string>> Metadata { get; set; }

        #endregion

        #region Public Properties

        public string ID { get; set; }

        public string Strand { get; set; }

        public Gene Gene { get; set; }

        public List<Exon> Exons { get; set; } = new List<Exon>();

        public List<Exon> CodingDomainSequences { get; set; } = new List<Exon>();

        public string ProteinID { get; set; }

        #endregion Public Properties

        #region Public Constructors

        public Transcript(string ID, Gene gene, MetadataListItem<List<string>> metadata, string ProteinID = null)
        {
            this.ID = ID;
            this.ProteinID = ProteinID == null ? ID : ProteinID;
            this.Gene = gene;
            this.Metadata = metadata;
            this.Strand = metadata.SubItems["strand"][0];
        }

        #endregion Public Constructors

        #region Public Methods

        public static List<string> combinatoricFailures = new List<string>();
        public IEnumerable<Protein> Translate(bool translateCodingDomains, bool includeVariants)
        {
            List<TranscriptPossiblyWithVariants> transcriptHaplotypes = CombineExonSequences(translateCodingDomains, includeVariants, out bool successfulCombination).Where(t => t.OkayToTranslate()).ToList();
            if (!successfulCombination)
            {
                combinatoricFailures.Add(ID);
                Console.WriteLine("combining exons failed " + combinatoricFailures.Count.ToString());
            }
            return ProteinAnnotation.OneFrameTranslationWithAnnotation(transcriptHaplotypes);
        }

        public IEnumerable<Protein> TranslateUsingAnnotatedStartCodons(Dictionary<Tuple<string, string, long>, List<Exon>> binnedCodingStarts, int indexBinSize, int minLength, bool includeVariants)
        {
            List<Exon> annotatedStarts = new List<Exon>();
            for (long i = Exons.Min(x => x.OneBasedStart) / indexBinSize; i < Exons.Max(x => x.OneBasedEnd) + 1; i++)
            {
                if (binnedCodingStarts.TryGetValue(new Tuple<string, string, long>(Gene.ChromID, Strand, i * indexBinSize), out List<Exon> exons))
                    annotatedStarts.AddRange(exons.Where(x => 
                        Exons.Any(xx => xx.Includes(Strand == "+" ? x.OneBasedStart : x.OneBasedEnd) // must include the start of the stop codon
                            && xx.Includes(Strand == "+" ? x.OneBasedStart + 2 : x.OneBasedEnd - 2)))); // and the end of the stop codon
            }

            char terminatingCharacter = ProteinAlphabet.Instance.GetFriendlyName(Alphabets.Protein.Ter)[0];
            if (annotatedStarts.Count > 0)
            {
                // gets the first annotated start that produces
                Dictionary<string, Protein> proteinDictionary = new Dictionary<string, Protein>();
                foreach (Exon annotatedStart in annotatedStarts)
                {
                    long startCodonStart = Strand == "+" ? annotatedStart.OneBasedStart : annotatedStart.OneBasedEnd; // CDS on the reverse strand have start and end switched
                    List<TranscriptPossiblyWithVariants> transcripts = CombineExonSequences(false, includeVariants, out bool success).Where(t => t.OkayToTranslate()).ToList();
                    foreach (TranscriptPossiblyWithVariants transcript in transcripts)
                    {
                        if (Strand == "+")
                        {
                            long exonLengthBeforeCodingStart = Exons.Where(x => x.OneBasedEnd < annotatedStart.OneBasedStart).Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                            long exonZeroBasedCodingStart = startCodonStart - Exons.FirstOrDefault(x => x.Includes(annotatedStart.OneBasedStart)).OneBasedStart;
                            transcript.ZeroBasedCodingStart = exonLengthBeforeCodingStart + exonZeroBasedCodingStart;
                            long lengthAfterCodingStart = transcript.Sequence.Count - transcript.ZeroBasedCodingStart;
                            transcript.Sequence = transcript.Sequence.GetSubSequence(transcript.ZeroBasedCodingStart, lengthAfterCodingStart);
                        }
                        else
                        {
                            long length = Exons.Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                            long chop = Exons.Where(x => x.OneBasedEnd >= annotatedStart.OneBasedEnd).Sum(x => annotatedStart.OneBasedEnd < x.OneBasedStart ? x.OneBasedEnd - x.OneBasedStart + 1 : x.OneBasedEnd - annotatedStart.OneBasedEnd);
                            long lengthAfterCodingStart = length - chop;
                            transcript.Sequence = transcript.Sequence.GetSubSequence(0, lengthAfterCodingStart);
                        }
                        Protein p = ProteinAnnotation.OneFrameTranslationWithAnnotation(transcript);
                        if (p.BaseSequence.Length >= minLength && ProteinAnnotation.AddProteinIfLessComplex(proteinDictionary, p))
                            yield return p;
                    }
                }
            }
            //return Translation.ThreeFrameTranslation(Exons, ProteinID);
        }

        /// <summary>
        /// Gets the possible exon haplotypes and combines them in all possible ways.
        /// </summary>
        /// <param name="exons"></param>
        /// <param name="includeVariants"></param>
        /// <returns></returns>
        public IEnumerable<TranscriptPossiblyWithVariants> CombineExonSequences(bool translateCodingDomains, bool includeVariants, out bool success, int maxCombosPerTranscript = 32)
        {
            int maxCombos = (int)Math.Log(maxCombosPerTranscript, 2) + 1;
            List<Exon> exons = translateCodingDomains ? CodingDomainSequences : Exons;

            // keep pruning until we get a reasonable amount of branching
            List<List<Exon>> exonSequences = null;
            int totalBranches = -1;
            int maxCombosForExons = maxCombos;
            while (totalBranches < 0 || totalBranches > maxCombosPerTranscript)
            {
                if (exons.Any(x => x.Variants.Sum(v => Variant.ParseVariantContext(v).Count(vv => vv.AlleleFrequency < 0.9)) > maxCombos))
                {
                    success = false;
                    return new List<TranscriptPossiblyWithVariants>();
                }
                exonSequences = exons.Select(x => x.GetExonSequences(maxCombosForExons, includeVariants, 0.9, false, 101)).ToList();
                totalBranches = (int)Math.Pow(2, exonSequences.Sum(possibleExons => possibleExons.Count - 1)) - Convert.ToInt32(exonSequences.Count == 0);
                maxCombosForExons = (int)Math.Ceiling((double)maxCombosForExons / (double)2);
                if (maxCombosForExons == 1 && totalBranches > maxCombosPerTranscript)
                {
                    // too hard of a problem for now without long-read sequencing data
                    success = false;
                    return new List<TranscriptPossiblyWithVariants>();
                }
            }
            totalBranches = (int)Math.Pow(2, exonSequences.Sum(possibleExons => possibleExons.Count - 1)) - Convert.ToInt32(exonSequences.Count == 0);
            TranscriptPossiblyWithVariants[] haplotypicSequences = new TranscriptPossiblyWithVariants[totalBranches];
            byte nByte = Encoding.UTF8.GetBytes("N".ToArray()).First();
            foreach (List<Exon> branchedExons in exonSequences)
            {
                for (int branch = 0; branch < branchedExons.Count; branch++)
                {
                    // using byte arrays for intermediate data structures greatly improves performance because of copying and storage efficiency over strings
                    byte[] branchedExonBytes = new byte[branchedExons[branch].Sequence.Count];
                    (branchedExons[branch].Sequence as Sequence).CopyTo(branchedExonBytes, 0, branchedExonBytes.Length);

                    int partitionLength = haplotypicSequences.Length / branchedExons.Count;
                    int startIdx = branch * partitionLength;
                    int endIdx = startIdx + partitionLength;
                    for (int j = startIdx; j < endIdx; j++)
                    {
                        if (haplotypicSequences[j] == null)
                        {
                            haplotypicSequences[j] = new TranscriptPossiblyWithVariants(this,
                                translateCodingDomains,
                                branchedExons[branch].Sequence,
                                branchedExons[branch].Sequence.Contains(nByte),
                                branchedExons[branch].Variants.OfType<Variant>().ToList());
                        }
                        else
                        {
                            byte[] newSeq = new byte[haplotypicSequences[j].Sequence.Count + branchedExonBytes.Length];
                            (haplotypicSequences[j].Sequence as Sequence).CopyTo(newSeq, 0, haplotypicSequences[j].Sequence.Count);
                            Array.Copy(newSeq, haplotypicSequences[j].Sequence.Count, branchedExonBytes, 0, branchedExonBytes.Length);
                        }
                    }
                }
            }
            success = true;
            return haplotypicSequences;
        }

        #endregion Public Methods

    }
}
