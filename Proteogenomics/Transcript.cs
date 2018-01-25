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

        /// <summary>
        /// This comes straight from the Bio.NET object
        /// </summary>
        private MetadataListItem<List<string>> Metadata { get; set; }

        #endregion

        #region Public Properties

        /// <summary>
        /// The transcript ID
        /// </summary>
        public string ID { get; set; }

        /// <summary>
        /// The strand this transcript lies on. Either '+' or '-'; might also be '.' if it's unknown, but would need to check.
        /// </summary>
        public string Strand { get; set; }

        /// <summary>
        /// The parent gene containing this transcript.
        /// </summary>
        public Gene Gene { get; set; }

        /// <summary>
        /// A list of exons, possibly including untranslated regions (UTRs).
        /// </summary>
        public List<Exon> Exons { get; set; } = new List<Exon>();

        /// <summary>
        /// A list of coding domain sequences (CDS), which are represented as Exons.
        /// </summary>
        public List<Exon> CodingDomainSequences { get; set; } = new List<Exon>();

        /// <summary>
        /// Variants annotated for this transcript using SnpEff.
        /// </summary>
        public List<SnpEffAnnotation> SnpEffVariants { get; set; } = new List<SnpEffAnnotation>();

        /// <summary>
        /// The protein ID derived from coding transcripts; this is imported from the GFF3 file if available.
        /// </summary>
        public string ProteinID { get; set; }

        #endregion Public Properties

        #region Translation Information

        public static List<string> combinatoricFailures = new List<string>();

        #endregion Translation Information

        #region Public Constructors

        /// <summary>
        /// Constructor from the GFF3 reader information, including IDs, strand and Protein ID if available.
        /// </summary>
        /// <param name="ID"></param>
        /// <param name="gene"></param>
        /// <param name="metadata"></param>
        /// <param name="ProteinID"></param>
        public Transcript(string ID, Gene gene, MetadataListItem<List<string>> metadata, string ProteinID = null)
        {
            this.ID = ID;
            this.ProteinID = ProteinID == null ? ID : ProteinID;
            this.Gene = gene;
            this.Metadata = metadata;
            this.Strand = metadata.SubItems["strand"][0];
        }

        #endregion Public Constructors

        #region Translate and Replace Methods

        /// <summary>
        /// Stores accessions for checking that they are unique.
        /// </summary>
        static HashSet<string> accessions = new HashSet<string>();

        /// <summary>
        /// Translates a transcript into full-length variant protein sequences using SnpEff annotations to modify protein sequences.
        /// </summary>
        /// <param name="translateCodingDomains"></param>
        /// <param name="includeVariants"></param>
        /// <param name="badProteinAccessions"></param>
        /// <param name="selenocysteineContaining"></param>
        /// <returns></returns>
        public IEnumerable<Protein> TranslateFromSnpEffAnnotatedSNVs(bool translateCodingDomains, bool includeVariants, HashSet<string> badProteinAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            badProteinAccessions = badProteinAccessions != null ? badProteinAccessions : new HashSet<string>();
            selenocysteineContaining = selenocysteineContaining != null ? selenocysteineContaining : new Dictionary<string, string>();
            Dictionary<string, Protein> proteinDictionary = new Dictionary<string, Protein>();

            // don't process proteins that have CDS without discrete sequence or without actual start or stop, i.e. containing 'X' or '*'
            if (badProteinAccessions.Contains(ProteinID))
                return new List<Protein>(); 

            Protein baseProtein = Translate(true, false, badProteinAccessions, selenocysteineContaining).FirstOrDefault();
            string proteinSequence = "";
            if (baseProtein != null)
                proteinSequence = baseProtein.BaseSequence;
            else
                return new List<Protein>();

            // todo: allow depth filter here using "DP" column, especially for indels
            double homozygousThreshold = 0.9;
            List<SnpEffAnnotation> synonymous = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> homozygousMissense = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> heterozygousMissense = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> other = new List<SnpEffAnnotation>();
            foreach (SnpEffAnnotation a in SnpEffVariants.Where(annotation => annotation.FeatureID == ID))
            {
                if (a.Synonymous) synonymous.Add(a);
                else if (a.Missense && a.Variant.Variants.FirstOrDefault(v => v.AlternateAllele == a.Allele).AlleleFrequency >= homozygousThreshold) homozygousMissense.Add(a);
                else if (a.Missense && a.Variant.Variants.FirstOrDefault(v => v.AlternateAllele == a.Allele).AlleleFrequency < homozygousThreshold) heterozygousMissense.Add(a);
                else other.Add(a);
            }

            int maxCombosPerTranscript = 32;
            List<List<SnpEffAnnotation>> missenseCombinations = new List<List<SnpEffAnnotation>> { new List<SnpEffAnnotation>(homozygousMissense) };
            missenseCombinations.AddRange(
                Enumerable.Range(1, heterozygousMissense.Count).SelectMany(k =>
                    ExtensionMethods.Combinations(heterozygousMissense, k)
                        .Select(heteroAlleles => new List<SnpEffAnnotation>(homozygousMissense).Concat(heteroAlleles).ToList()))
                    .ToList());

            // fake phasing (of variants within 100 bp of each other) to eliminate combinitorial explosion
            // this section eliminates instances where variants that are very close to each other are emitted as separate haplotypes
            int fakePhasingRange = 101;
            int range = fakePhasingRange;
            List<List<SnpEffAnnotation>> fakePhased = null;
            while (fakePhased == null || fakePhased.Count > Math.Max(2, maxCombosPerTranscript))
            {
                fakePhased = new List<List<SnpEffAnnotation>>();
                List<List<SnpEffAnnotation>> phasedLast = missenseCombinations.OrderBy(x => x.Count).ToList();
                List<int> variantSites = phasedLast.Last().Select(v => v.Variant.Start).ToList();
                foreach (List<SnpEffAnnotation> hap in phasedLast)
                {
                    List<int> hapSites = hap.Select(v => v.Variant.Start).ToList();
                    bool containsVariantsToFakePhase = variantSites.Any(vs => !hapSites.Contains(vs) && hap.Any(v => Math.Abs(v.Variant.Start - vs) < range));
                    if (!containsVariantsToFakePhase)
                        fakePhased.Add(hap);
                }
                missenseCombinations = fakePhased;
                range *= 2;
            }

            List<Protein> proteins = new List<Protein>();
            foreach (List<SnpEffAnnotation> annotations in missenseCombinations)
            {
                string variantAminoAcidSequence = proteinSequence;
                foreach (SnpEffAnnotation a in annotations)
                {
                    variantAminoAcidSequence = variantAminoAcidSequence.Substring(0, a.AminoAcidLocation - 1) + a.AlternateAminoAcid + variantAminoAcidSequence.Substring(a.AminoAcidLocation, variantAminoAcidSequence.Length - a.AminoAcidLocation);
                }
                List<SnpEffAnnotation> combinedAnnotations = synonymous.Concat(annotations).ToList();
                string proteinSnpEffAnnotation = "{" + String.Join(" ", combinedAnnotations.Select(a => "AF=" + a.Variant.Variants.FirstOrDefault(v => v.AlternateAllele == a.Allele).AlleleFrequency.ToString("N2") + ";" + a.Annotation)) + "} OS=Homo sapiens GN=" + Gene.ID;
                List<SequenceVariation> sequenceVariations = combinedAnnotations.Select(a => new SequenceVariation(a.Variant.Start, a.Variant.Reference.BaseString, a.Allele, "AF=" + a.Variant.Variants.FirstOrDefault(v => v.AlternateAllele == a.Allele).AlleleFrequency.ToString("N2") + ";ANN=" + a.Annotation)).ToList();
                string accession = ProteinID;
                int arbitraryNumber = 1;
                while (accessions.Contains(accession))
                {
                    accession = ProteinID + "_v" + arbitraryNumber++.ToString();
                }
               proteins.Add(new Protein(variantAminoAcidSequence.Split('*')[0], accession, null, null, null, null, proteinSnpEffAnnotation, false, false, null, sequenceVariations));
            }
            return proteins;
        }

        #endregion Translate and Replace Methods

        #region Translate from Variant Nucleotide Sequences Methods

        public IEnumerable<Protein> Translate(bool translateCodingDomains, bool includeVariants, HashSet<string> badProteinAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            List<TranscriptPossiblyWithVariants> transcriptHaplotypes = CombineExonSequences(translateCodingDomains, includeVariants, out bool successfulCombination).Where(t => t.OkayToTranslate()).ToList();
            if (!successfulCombination)
            {
                combinatoricFailures.Add(ID);
                Console.WriteLine("combining exons failed " + combinatoricFailures.Count.ToString());
            }
            return ProteinAnnotation.OneFrameTranslationWithAnnotation(transcriptHaplotypes, badProteinAccessions, selenocysteineContaining);
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
                            long lengthAfterCodingStart = transcript.VariantTranscriptSequence.Count - transcript.ZeroBasedCodingStart;
                            transcript.VariantTranscriptSequence = transcript.VariantTranscriptSequence.GetSubSequence(transcript.ZeroBasedCodingStart, lengthAfterCodingStart);
                        }
                        else
                        {
                            long length = Exons.Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                            long chop = Exons.Where(x => x.OneBasedEnd >= annotatedStart.OneBasedEnd).Sum(x => annotatedStart.OneBasedEnd < x.OneBasedStart ? x.OneBasedEnd - x.OneBasedStart + 1 : x.OneBasedEnd - annotatedStart.OneBasedEnd);
                            long lengthAfterCodingStart = length - chop;
                            transcript.VariantTranscriptSequence = transcript.VariantTranscriptSequence.GetSubSequence(0, lengthAfterCodingStart);
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
            if (Strand != "+")
                exonSequences.Reverse();
            foreach (List<Exon> branchedExons in exonSequences)
            {
                for (int branch = 0; branch < branchedExons.Count; branch++)
                {
                    // using byte arrays for intermediate data structures greatly improves performance because of copying and storage efficiency over strings
                    Sequence sequence = Strand == "+" ? branchedExons[branch].Sequence as Sequence : branchedExons[branch].Sequence.GetReverseComplementedSequence() as Sequence;
                    byte[] branchedExonBytes = new byte[branchedExons[branch].Sequence.Count];
                    sequence.CopyTo(branchedExonBytes, 0, branchedExonBytes.Length);
                    string seq = new string(branchedExonBytes.Select(b => char.ToUpperInvariant((char)b)).ToArray()); // debugging

                    int partitionLength = haplotypicSequences.Length / branchedExons.Count;
                    int startIdx = branch * partitionLength;
                    int endIdx = startIdx + partitionLength;
                    for (int j = startIdx; j < endIdx; j++)
                    {
                        if (haplotypicSequences[j] == null)
                        {
                            haplotypicSequences[j] = new TranscriptPossiblyWithVariants(this,
                                translateCodingDomains,
                                sequence,
                                branchedExons[branch].Sequence.Contains(nByte),
                                branchedExons[branch].Variants.OfType<Variant>().ToList());
                        }
                        else
                        {
                            byte[] newSeq = new byte[haplotypicSequences[j].VariantTranscriptSequence.Count + branchedExonBytes.Length];
                            (haplotypicSequences[j].VariantTranscriptSequence as Sequence).CopyTo(newSeq, 0, haplotypicSequences[j].VariantTranscriptSequence.Count);
                            Array.Copy(branchedExonBytes, 0, newSeq, haplotypicSequences[j].VariantTranscriptSequence.Count, branchedExonBytes.Length);
                            haplotypicSequences[j].VariantTranscriptSequence = new Sequence(haplotypicSequences[j].VariantTranscriptSequence.Alphabet, newSeq);
                        }
                    }
                }
            }
            success = true;
            return haplotypicSequences;
        }

        #endregion Translate from Variant Nucleotide Sequences Methods

    }
}
