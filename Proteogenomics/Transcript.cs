using Bio;
using Bio.Extensions;
using Bio.VCF;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class Transcript
        : Interval
    {
        #region Public Properties

        /// <summary>
        /// The transcript ID
        /// </summary>
        public string ID { get; set; }

        /// <summary>
        /// The transcript version
        /// </summary>
        public string Version { get; set; }

        /// <summary>
        /// The parent gene containing this transcript.
        /// </summary>
        public Gene Gene { get; set; }

        /// <summary>
        /// A list of exons, possibly including untranslated regions (UTRs).
        /// </summary>
        public List<Exon> Exons { get; set; } = new List<Exon>();

        /// <summary>
        /// Exons sorted by start position or by reverse end position if on reverse strand
        /// </summary>
        public List<Exon> ExonsSortedStrand { get { return SortedStrand(); } }

        /// <summary>
        /// List of the untranslated regions
        /// </summary>
        public List<UTR> UTRs { get; set; } = new List<UTR>();

        /// <summary>
        /// Downstream region, 5 kb by default
        /// </summary>
        public Downstream Downstream { get; set; }

        /// <summary>
        /// Upstream region, 5 kb by default
        /// </summary>
        public Upstream Upstream { get; set; }

        /// <summary>
        /// Introns for this transcript
        /// </summary>
        public List<Intron> Introns { get; set; }

        /// <summary>
        /// Coding domain sequence (CDS) information
        /// </summary>
        public List<CDS> CodingDomainSequences { get; set; } = new List<CDS>();

        /// <summary>
        /// Coding sequence
        /// </summary>
        public ISequence CodingSequence { get; set; }

        /// <summary>
        /// Variants annotated for this transcript using SnpEff.
        /// </summary>
        public List<SnpEffAnnotation> SnpEffVariants { get; set; } = new List<SnpEffAnnotation>();

        /// <summary>
        /// The protein ID derived from coding transcripts; this is imported from the GFF3 file if available.
        /// </summary>
        public string ProteinID { get; set; }

        /// <summary>
        /// Annotations for the variants applied to this transcript
        /// </summary>
        public List<string> VariantAnnotations { get; set; } = new List<string>();

        #endregion Public Properties

        #region Translation Information

        public const int DEFAULT_UP_DOWN_LENGTH = 5000;

        public static List<string> combinatoricFailures = new List<string>();

        #endregion Translation Information

        #region Public Constructors

        /// <summary>
        /// Constructor from the GFF3 reader information, including IDs, strand and Protein ID if available.
        /// </summary>
        /// <param name="id"></param>
        /// <param name="gene"></param>
        /// <param name="metadata"></param>
        /// <param name="ProteinID"></param>
        public Transcript(string id, string version, Gene gene, string strand, long oneBasedStart, long oneBasedEnd, string ProteinID = null)
            : base(gene.ChromosomeID, strand, oneBasedStart, oneBasedEnd)
        {
            this.ID = id;
            this.Version = version;
            this.ProteinID = ProteinID ?? id;
            this.Gene = gene;
        }

        #endregion Public Constructors

        #region Translate and Replace with SnpEff Annotated Variations Method

        /// <summary>
        /// Apply this variant and adjust the start and stop indices
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        public override Interval ApplyVariant(Variant variant)
        {
            Interval interval = base.ApplyVariant(variant);
            Transcript transcript = new Transcript(ID, Version, Gene, interval.Strand, interval.OneBasedStart, interval.OneBasedEnd, ProteinID);
            for (int i = 0; i < CodingDomainSequences.Count; i++)
            {
                if (CodingDomainSequences[i].Includes(variant))
                {
                    transcript.CodingDomainSequences.Add(CodingDomainSequences[i].ApplyVariant(variant) as CDS);
                }
                transcript.CodingDomainSequences.Add(CodingDomainSequences[i]);
            }
            for (int i = 0; i < Exons.Count; i++)
            {
                if (Exons[i].Includes(variant))
                {
                    transcript.Exons.Add(Exons[i].ApplyVariant(variant) as Exon); // applies variant to sequence
                }
                transcript.Exons.Add(Exons[i]);
            }

            transcript.VariantAnnotations = new List<string>(VariantAnnotations);
            transcript.CreateIntrons();
            transcript.CreateUTRs();
            transcript.CreateUpDown(Gene.Chromosome);
            return transcript;
        }

        /// <summary>
        /// Takes in list of variants to combinitorially apply to this transcript
        /// </summary>
        /// <param name="variantOrderedDescStart"></param>
        /// <returns></returns>
        public IEnumerable<Transcript> ApplyVariantsCombinitorially(List<Variant> variantOrderedDescStart)
        {
            List<Transcript> result = new List<Transcript>();
            for (int i = 0; i < variantOrderedDescStart.Count; i++)
            {
                Variant v = variantOrderedDescStart[i];
                VariantEffect variantEffect = AnnotateVariant(v);
                if (variantEffect == null)
                {
                    continue;
                }
                Transcript newTranscript = ApplyVariant(v) as Transcript;
                newTranscript.VariantAnnotations.Add(variantEffect.ToString());
                result.Add(newTranscript);
                if (variantOrderedDescStart.Count - i > 1)
                {
                    result.AddRange(newTranscript.ApplyVariantsCombinitorially(variantOrderedDescStart.GetRange(i + 1, variantOrderedDescStart.Count - i - 1)));
                }
                if (v.GenotypeType == GenotypeType.HETEROZYGOUS)
                {
                    result.Add(this);
                }
                if (v.GenotypeType == GenotypeType.HETEROZYGOUS && variantOrderedDescStart.Count - i > 1)
                {
                    result.AddRange(this.ApplyVariantsCombinitorially(variantOrderedDescStart.GetRange(i + 1, variantOrderedDescStart.Count - i - 1)));
                }
            }
            return result;
        }

        /// <summary>
        /// Gets a string representing a variant applied to this transcript 
        /// </summary>
        /// <param name="variant"></param>
        public VariantEffect AnnotateVariant(Variant variant)
        {
            // Leaving off here for 180319


            // >> Figure out the codon change from the information in this transcript (port from SnpEff CodonChange)
            // >> Figure out the VariantEffect for the change, and note the high/medium/low/modifier status of the change
            // Translate the codons like I was doing in ProteinAnnotation, making that byte array to stick into translate. There are places in CodonChange and VariantEffect that need this

            //return "variant:"; // frameshift or something SPACE notations about what changes were made to the sequence

            // Then in translation, make a simple method to look up the bad accessions
            //      (could also try to assess this from the sequence itself using the warnings and errors)
            //      (namely, does it have stop codons in it, does it have a start codon)
            //      (but will have to do this anyway to find the selenocysteine sequences, so might as well just keep that code)
            // Put together coding sequences based on the UTRs, now that they're being fixed in the ApplyVariants method
            // Translate the thing
            //      Keep the accessions to look up and increment accessions to make them unique during translation

            // Delete all the ProteinAnnotation code and Transcript and TranscriptPossiblyWithVariants code that isn't being used
            //     Don't need to prepare for translation anymore
            //     Don't need to step through translation anymore (but keep the formatting of variants from ProteinAnnotation)
            //  - space delimited with variant annotations: "variant:type originalSequence###alteredSequence chr:oneBasedStart"
            //     Don't need to TranslateFromSnpEff 

            // Test
            //   1. UTR ranges get change
            //   2. Correct UTR gets changed (5' or 3')
            //   3. Variants get applied correctly
            //   Uh, lots more.


            if (!Intersects(variant)) return null; // Sanity check

            // Large structural variant including the whole transcript?
            if (variant.Includes(this) && variant.isStructural())
            {
                CodonChange codonChange = CodonChange.Factory(variant, this, variantEffects);
                codonChange.ChangeCodon();
                return true;
            }

            //---
            // Structural variants may affect more than one exon
            //---
            bool mayAffectSeveralExons = variant.isStructural() || variant.isMixed() || variant.isMnp();
            if (mayAffectSeveralExons)
            {
                int countExon = 0;
                foreach (Exon ex in Exons)
                {
                    if (ex.Intersects(variant))
                    {
                        countExon++;
                    }
                }

                // More than one exon?
                if (countExon > 1)
                {
                    CodonChange codonChange = CodonChange.Factory(variant, this, variantEffects);
                    codonChange.ChangeCodon();
                    return true;
                }
            }

            //---
            // Does it hit an exon?
            // Note: This only adds spliceSites effects, for detailed codon
            //       changes effects we use 'CodonChange' class
            //---
            bool exonAnnotated = false;
            foreach (Exon ex in Exons)
            {
                if (ex.Intersects(variant))
                {
                    exonAnnotated |= ex.variantEffect(variant, variantEffects);
                }
            }

            //---
            // Hits a UTR region?
            //---
            bool included = false;
            for (Utr utr : utrs)
            {
                if (utr.intersects(variant))
                {
                    // Calculate the effect
                    utr.variantEffect(variant, variantEffects);
                    included |= utr.includes(variant); // Is this variant fully included in the UTR?
                }
            }
            if (included)
            {
                return true; // Variant fully included in the UTR? => We are done.
            }

            //---
            // Does it hit an intron?
            //---
            foreach (Intron intron in Introns)
            {
                if (intron.Intersects(variant))
                {
                    intron.variantEffect(variant, variantEffects);
                    included |= intron.Includes(variant); // Is this variant fully included in this intron?
                }
            }
            if (included)
            {
                return true; // Variant fully included? => We are done.
            }

            //---
            // No annotations from exons? => Add transcript
            //---
            if (!exonAnnotated)
            {
                variantEffects.Add(variant, this, EffectType.TRANSCRIPT, "");
                return true;
            }

            return exonAnnotated;
        }

        /// <summary>
        /// Stores accessions for checking that they are unique.
        /// </summary>
        private static HashSet<string> accessions = new HashSet<string>();

        /// <summary>
        /// Translates a transcript into full-length variant protein sequences using SnpEff annotations to modify protein sequences.
        /// </summary>
        /// <param name="translateCodingDomains"></param>
        /// <param name="includeVariants"></param>
        /// <param name="badProteinAccessions"></param>
        /// <param name="selenocysteineContaining"></param>
        /// <returns></returns>
        public IEnumerable<Protein> TranslateFromSnpEffAnnotatedSNVs(bool translateCodingDomains, bool includeVariants, string reference, Dictionary<string, string> proteinSequences, HashSet<string> badProteinAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            badProteinAccessions = badProteinAccessions != null ? badProteinAccessions : new HashSet<string>();
            selenocysteineContaining = selenocysteineContaining != null ? selenocysteineContaining : new Dictionary<string, string>();

            // don't process proteins that have CDS without discrete sequence or without actual start or stop, i.e. containing 'X' or '*'
            proteinSequences.TryGetValue(ProteinID, out string proteinSequence);
            if (proteinSequence == null || badProteinAccessions.Contains(ProteinID))
            {
                return new List<Protein>();
            }

            // todo: allow depth filter here using "DP" column, especially for indels
            List<SnpEffAnnotation> synonymous = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> homozygousMissense = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> heterozygousMissense = new List<SnpEffAnnotation>();
            List<SnpEffAnnotation> other = new List<SnpEffAnnotation>();
            string id = reference.StartsWith("GRCh38") ? ID + "." + Version : ID;
            foreach (SnpEffAnnotation a in SnpEffVariants.Where(annotation => annotation.FeatureID == id))
            {
                if (a.Synonymous)
                {
                    synonymous.Add(a);
                }
                else if (a.Missense && a.Variant.Variants.FirstOrDefault(v => v.AlternateAlleleString == a.Allele).GenotypeType == GenotypeType.HOMOZYGOUS_ALT)
                {
                    homozygousMissense.Add(a);
                }
                else if (a.Missense && a.Variant.Variants.FirstOrDefault(v => v.AlternateAlleleString == a.Allele).GenotypeType == GenotypeType.HETEROZYGOUS)
                {
                    heterozygousMissense.Add(a);
                }
                else
                {
                    other.Add(a);
                }
            }

            int maxCombosPerTranscript = 32;
            List<List<SnpEffAnnotation>> missenseCombinations = new List<List<SnpEffAnnotation>> { new List<SnpEffAnnotation>(homozygousMissense) };
            missenseCombinations.AddRange(
                Enumerable.Range(1, heterozygousMissense.Count).SelectMany(k =>
                    ProteogenomicsUtility.Combinations(heterozygousMissense, k)
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
                List<long> variantSites = phasedLast.Last().Select(v => v.Variant.OneBasedStart).ToList();
                foreach (List<SnpEffAnnotation> hap in phasedLast)
                {
                    List<long> hapSites = hap.Select(v => v.Variant.OneBasedStart).ToList();
                    bool containsVariantsToFakePhase = variantSites.Any(vs => !hapSites.Contains(vs) && hap.Any(v => Math.Abs(v.Variant.OneBasedStart - vs) < range));
                    if (!containsVariantsToFakePhase)
                    {
                        fakePhased.Add(hap);
                    }
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
                List<SequenceVariation> sequenceVariations = combinedAnnotations.Select(a => new SequenceVariation(a.Variant.OneBasedStart, a.Variant.ReferenceAlleleString, a.Allele, "AF=" + a.Variant.Variants.FirstOrDefault(v => v.AlternateAllele == a.Allele).AlleleFrequency.ToString("N2") + ";ANN=" + a.Annotation)).ToList();
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

        #endregion Translate and Replace with SnpEff Annotated Variations Method

        #region Translate from Variant Nucleotide Sequences Methods

        /// <summary>
        /// Get the coding sequence for this transcript.
        /// SnpEff keeps track of the UTRs to figure this out. I suppose that will work, now that I'm using the interval tree to dive down to change those ranges.
        /// </summary>
        /// <returns></returns>
        public ISequence RetrieveCodingSequence()
        {
            if (CodingSequence != null)
            {
                return CodingSequence;
            }

            // Concatenate all exons
            List<Exon> exons = ExonsSortedStrand;
            StringBuilder sequence = new StringBuilder();
            int utr5len = 0, utr3len = 0;

            // 5 prime UTR length
            foreach (UTR utr in UTRs.OfType<UTR5Prime>())
            {
                utr5len += (int)utr.Length();
            }

            // Append all exon sequences
            IAlphabet alphabet = null;
            bool missingSequence = false;
            foreach (Exon exon in exons)
            {
                missingSequence |= exon.Sequence != null; // If there is no sequence, we are in trouble
                sequence.Append(Strand == "+" ? exon.Sequence : exon.Sequence.GetReverseComplementedSequence());
                alphabet = alphabet == null || alphabet.HasAmbiguity && !exon.Sequence.Alphabet.HasAmbiguity ? // keep the alphabet with the most characters
                    alphabet :
                    exon.Sequence.Alphabet;
            }

            if (missingSequence)
            {
                CodingSequence = new Sequence(Alphabets.DNA, ""); // One or more exons does not have sequence. Nothing to do
            }
            else
            {
                // OK, all exons have sequences

                // 3 prime UTR length
                foreach (UTR utr in UTRs.OfType<UTR3Prime>())
                {
                    utr3len += (int)utr.Length();
                }

                // Cut 5 prime UTR and 3 prime UTR points
                int subEnd = sequence.Length - utr3len;

                if (utr5len > subEnd)
                {
                    CodingSequence = new Sequence(Alphabets.DNA, "");
                }
                else
                {
                    CodingSequence = new Sequence(alphabet, sequence.ToString().Substring(utr5len, subEnd));
                }
            }
            return CodingSequence;
        }

        /// <summary>
        /// Get strands sorted by start (forward) or end and in reverse (reverse strand)
        /// </summary>
        /// <returns></returns>
        private List<Exon> SortedStrand()
        {
            return ExonsSortedStrand != null ?
                ExonsSortedStrand :
                Strand == "+" ?
                    Exons.OrderBy(x => x.OneBasedStart).ToList() :
                    Exons.OrderByDescending(x => x.OneBasedEnd).ToList();
        }

        public IEnumerable<Protein> Translate(bool translateCodingDomains, bool includeVariants, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            List<TranscriptPossiblyWithVariants> transcriptHaplotypes = CombineExonSequences(translateCodingDomains, includeVariants, out bool successfulCombination).Where(t => t.OkayToTranslate()).ToList();
            if (!successfulCombination)
            {
                combinatoricFailures.Add(ID);
                Console.WriteLine("combining exons failed " + combinatoricFailures.Count.ToString());
            }
            return ProteinAnnotation.OneFrameTranslationWithAnnotation(transcriptHaplotypes, incompleteTranscriptAccessions, selenocysteineContaining);
        }

        public IEnumerable<Protein> TranslateUsingAnnotatedStartCodons(Dictionary<Tuple<string, string, long>, List<CDS>> binnedCodingStarts, int indexBinSize, int minLength, bool includeVariants)
        {
            List<CDS> annotatedStarts = new List<CDS>();
            for (long i = Exons.Min(x => x.OneBasedStart) / indexBinSize; i < Exons.Max(x => x.OneBasedEnd) + 1; i++)
            {
                if (binnedCodingStarts.TryGetValue(new Tuple<string, string, long>(Gene.ChromosomeID, Strand, i * indexBinSize), out List<CDS> cds))
                {
                    annotatedStarts.AddRange(cds.Where(x =>
                        Exons.Any(xx => xx.Includes(Strand == "+" ? x.OneBasedStart : x.OneBasedEnd) // must include the start of the stop codon
                            && xx.Includes(Strand == "+" ? x.OneBasedStart + 2 : x.OneBasedEnd - 2)))); // and the end of the stop codon
                }
            }

            char terminatingCharacter = ProteinAlphabet.Instance.GetFriendlyName(Alphabets.Protein.Ter)[0];
            if (annotatedStarts.Count > 0)
            {
                // gets the first annotated start that produces
                Dictionary<string, Protein> proteinDictionary = new Dictionary<string, Protein>();
                foreach (CDS annotatedStart in annotatedStarts)
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
                        {
                            yield return p;
                        }
                    }
                }
            }
            //return Translation.ThreeFrameTranslation(Exons, ProteinID);
        }

        public bool isProteinCoding()
        {
            return CodingDomainSequences.Count > 0;
        }

        #endregion Translate from Variant Nucleotide Sequences Methods

        #region Create Interval Methods

        /// <summary>
        /// Create UTR regions for this transcript
        /// </summary>
        public List<UTR> CreateUTRs()
        {
            if (CodingDomainSequences.Count == 0)
                return UTRs;

            CodingDomainSequences = CodingDomainSequences.OrderBy(c => c.OneBasedStart).ToList();
            long codingLeft = CodingDomainSequences.First().OneBasedStart;
            long codingRight = CodingDomainSequences.Last().OneBasedEnd;
            foreach (Exon x in Exons.OrderBy(c => c.OneBasedStart))
            {
                if (x.OneBasedStart < codingLeft)
                {
                    long end = x.OneBasedEnd < codingLeft ? x.OneBasedEnd : codingLeft;
                    UTRs.Add(Strand == "+" ? new UTR5Prime(ChromosomeID, Strand, x.OneBasedStart, end) as UTR : new UTR3Prime(ChromosomeID, Strand, x.OneBasedStart, end) as UTR);
                }
                if (x.OneBasedEnd > CodingDomainSequences.Last().OneBasedEnd)
                {
                    long start = x.OneBasedStart < codingRight ? x.OneBasedStart : codingRight;
                    UTRs.Add(Strand == "+" ? new UTR3Prime(ChromosomeID, Strand, start, x.OneBasedEnd) as UTR : new UTR5Prime(ChromosomeID, Strand, start, x.OneBasedEnd) as UTR);
                }
            }
            return UTRs;
        }

        /// <summary>
        /// Creates upstream and downstream intervals for this transcript
        /// </summary>
        /// <param name="chromosomeSequence"></param>
        public List<Interval> CreateUpDown(Chromosome chromosomeSequence)
        {
            long chrMin = 1;
            long chrMax = chromosomeSequence.Sequence.Count;

            // Create up/down stream intervals and add them to the list
            long beforeStart = Math.Max(chrMin - DEFAULT_UP_DOWN_LENGTH, chrMin);
            long beforeEnd = Math.Max(chrMin - 1, chrMin);
            long afterStart = Math.Min(OneBasedEnd + 1, chrMax);
            long afterEnd = Math.Min(OneBasedEnd + DEFAULT_UP_DOWN_LENGTH, chrMax);

            if (Strand == "+")
            {
                if (beforeStart < beforeEnd) Upstream = new Upstream(chromosomeSequence.ChromosomeID, Strand, beforeStart, beforeEnd);
                if (afterStart < afterEnd) Downstream = new Downstream(chromosomeSequence.ChromosomeID, Strand, afterStart, afterEnd);
            }
            else
            {
                if (afterStart < afterEnd) Upstream = new Upstream(chromosomeSequence.ChromosomeID, Strand, afterStart, afterEnd);
                if (beforeStart < beforeEnd) Downstream = new Downstream(chromosomeSequence.ChromosomeID, Strand, beforeStart, beforeEnd);
            }

            return new List<Interval> { Upstream, Downstream };
        }

        public List<Intron> CreateIntrons()
        {
            Exon previous = null;
            foreach (Exon x in Exons)
            {
                if (previous == null)
                {
                    previous = x;
                    continue;
                }
                Intron intron = new Intron(x.ChromosomeID, x.Strand, previous.OneBasedEnd + 1, x.OneBasedStart - 1);
                if (intron.Length() > 0)
                    Introns.Add(intron);
            }
            return Introns;
        }

        #endregion Create Interval Methods

        #region Warning Methods

        /// <summary>
        /// Check if coding length is multiple of 3 in protein coding transcripts
        /// </summary>
        /// <returns></returns>
        public bool isErrorProteinLength(bool treatAllAsProteinCoding)
        {
            if (!treatAllAsProteinCoding && !isProteinCoding()) return false;
            return (cds().length() % 3) != 0;
        }

        /// <summary>
        /// Is the first codon a START codon?
        /// </summary>
        /// <returns></returns>
        public bool isErrorStartCodon()
        {
            if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return false;

            // Not even one codon in this protein? Error
            string cds = cds();
            if (cds.length() < 3) return true;

            string codon = cds.substring(0, 3);
            return !codonTable().isStart(codon);
        }

        /// <summary>
        /// Check if protein sequence has STOP codons in the middle of the coding sequence
        /// </summary>
        /// <returns></returns>
        public bool isErrorStopCodonsInCds()
        {
            //if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return false;

            // Get protein sequence
            string prot = protein();
            if (prot == null) return false;

            // Any STOP codon before the end?
            char bases[] = prot.toCharArray();
            int max = bases.length - 1;
            int countErrs = 0;
            for (int i = 0; i < max; i++)
                if (bases[i] == '*')
                {
                    countErrs++;
                    // We allow up to one STOP codon because it can be a RARE_AMINO_ACID which is coded as a STOP codon.
                    // More than one STOP codon is not "normal", so it's probably an error in the genomic annotations (e.g. ENSEMBL or UCSC)
                    if (countErrs > 1) return true;
                }

            // OK
            return false;
        }

        /// <summary>
        /// Is the last codon a STOP codon?
        /// </summary>
        /// <returns></returns>
        public bool isWarningStopCodon()
        {
            if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return false;

            // Not even one codon in this protein? Error
            string cds = cds();
            if (cds.length() < 3) return true;

            string codon = cds.substring(cds.length() - 3);
            return !codonTable().isStop(codon);
        }


        #endregion
    }
}