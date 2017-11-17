using Bio;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

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

        public IEnumerable<Protein> Translate(bool translateCodingDomains, bool includeVariants)
        {
            List<TranscriptPossiblyWithVariants> transcriptHaplotypes = CombineExonSequences(translateCodingDomains, includeVariants).Where(t => t.OkayToTranslate()).ToList();
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
                    List<TranscriptPossiblyWithVariants> transcripts = CombineExonSequences(false, includeVariants).Where(t => t.OkayToTranslate()).ToList();
                    foreach (TranscriptPossiblyWithVariants transcript in transcripts)
                    {
                        if (Strand == "+")
                        {
                            long exonLengthBeforeCodingStart = Exons.Where(x => x.OneBasedEnd < annotatedStart.OneBasedStart).Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                            long exonZeroBasedCodingStart = startCodonStart - Exons.FirstOrDefault(x => x.Includes(annotatedStart.OneBasedStart)).OneBasedStart;
                            transcript.ZeroBasedCodingStart = exonLengthBeforeCodingStart + exonZeroBasedCodingStart;
                            long lengthAfterCodingStart = transcript.Sequence.Length - transcript.ZeroBasedCodingStart;
                            transcript.Sequence = SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, transcript.Sequence).GetSubSequence(transcript.ZeroBasedCodingStart, lengthAfterCodingStart));
                        }
                        else
                        {
                            long length = Exons.Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                            long chop = Exons.Where(x => x.OneBasedEnd >= annotatedStart.OneBasedEnd).Sum(x => annotatedStart.OneBasedEnd < x.OneBasedStart ? x.OneBasedEnd - x.OneBasedStart + 1 : x.OneBasedEnd - annotatedStart.OneBasedEnd);
                            long lengthAfterCodingStart = length - chop;
                            transcript.Sequence = SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, transcript.Sequence).GetSubSequence(0, lengthAfterCodingStart));
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
        public List<TranscriptPossiblyWithVariants> CombineExonSequences(bool translateCodingDomains, bool includeVariants)
        {
            List<Exon> exons = translateCodingDomains ? CodingDomainSequences : Exons;
            List<List<Exon>> exonSequences = exons.Select(x => x.GetExonSequences(includeVariants, 0.9).ToList()).ToList();
            List<TranscriptPossiblyWithVariants> sequences = new List<TranscriptPossiblyWithVariants>();
            foreach (List<Exon> nextExon in exonSequences)
            {
                sequences = sequences.Count == 0 
                    ?
                    nextExon.Select(x => new TranscriptPossiblyWithVariants(this, translateCodingDomains, SequenceExtensions.ConvertToString(x.Sequence), x.Variants.OfType<Variant>().ToList())).ToList() 
                    :
                    (from curr in sequences
                    from next in nextExon
                    select new TranscriptPossiblyWithVariants(this, translateCodingDomains, curr.Sequence + SequenceExtensions.ConvertToString(next.Sequence), curr.Variants.Concat(next.Variants.OfType<Variant>()).ToList())
                        ).ToList();
            }
            return sequences;
        }

        #endregion Public Methods

    }
}
