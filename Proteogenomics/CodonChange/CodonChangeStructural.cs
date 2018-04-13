using Bio.Extensions;
using System;

namespace Proteogenomics
{
    public abstract class CodonChangeStructural
        : CodonChange
    {
        protected bool coding;
        protected int exonFull, exonPartial;
        protected string cdsAlt;
        protected string cdsRef;

        protected CodonChangeStructural(Variant variant, Transcript transcript, VariantEffects variantEffects)
            : base(variant, transcript, variantEffects)
        {
            coding = transcript.IsProteinCoding(); // || Config.get().isTreatAllAsProteinCoding();
            CountAffectedExons();
        }

        /// <summary>
        /// Differences between two CDSs after removing equal codons from
        /// the beginning and from the end of both strings
        /// </summary>
        protected void cdsDiff()
        {
            int min = Math.Min(cdsRef.Length, cdsAlt.Length) / 3;

            // Removing form the beginning
            CodonStartNumber = 0;
            CodonStartIndex = 0;
            for (int i = 0; i <= min; i++)
            {
                CodonStartNumber = i;

                // Codons differ?
                if (!codonEquals(cdsRef, cdsAlt, i, i))
                {
                    // Find index difference within codon
                    CodonStartIndex = codonDiffIndex(cdsRef, cdsAlt, i, i);
                    break;
                }
            }

            // Removing trailing codons
            int codonNumEndRef = cdsRef.Length / 3;
            int codonNumEndAlt = cdsAlt.Length / 3;

            for (; codonNumEndRef >= CodonStartNumber && codonNumEndAlt >= CodonStartNumber; codonNumEndRef--, codonNumEndAlt--)
            {
                if (!codonEquals(cdsRef, cdsAlt, codonNumEndRef, codonNumEndAlt))
                {
                    break;
                }
            }

            // Codons
            CodonsReference = Codons(cdsRef, CodonStartNumber, codonNumEndRef);
            CodonsAlternate = Codons(cdsAlt, CodonStartNumber, codonNumEndAlt);

            // No codon difference found?
            if (CodonsReference == "" && CodonsAlternate == "")
            {
                CodonStartNumber = -1;
                CodonStartIndex = -1;
            }
        }

        public override void ChangeCodon()
        {
            if (Variant.Includes(Transcript))
            {
                // Whole transcript affected?
                EffectTranscript();
            }
            else
            {
                // Does the variant affect any exons?
                if (exonFull > 0 || exonPartial > 0)
                {
                    Exons();
                }
                else
                {
                    Intron();
                }
            }
        }

        protected void codonChangeSuper()
        {
            base.ChangeCodon();
        }

        private static int codonDiffIndex(string cdsRef, string cdsAlt, int codonNumRef, int codonNumAlt)
        {
            for (int h = 0, i = 3 * codonNumRef, j = 3 * codonNumAlt; h < 3; i++, j++, h++)
            {
                // Premature end of sequence? (i.e. sequence ends before codon end)
                if (i >= cdsRef.Length || j >= cdsAlt.Length)
                {
                    return h;
                }

                // Different base? Return index
                if (cdsRef[i] != cdsAlt[j])
                {
                    return h;
                }
            }

            return -1;
        }

        /// <summary>
        /// Compare codons from cdsRef[codonNumRef] and cdsAlt[codonNumAlt]
        /// </summary>
        /// <param name="cdsRef"></param>
        /// <param name="cdsAlt"></param>
        /// <param name="codonNumRef"></param>
        /// <param name="codonNumAlt"></param>
        /// <returns></returns>
        private static bool codonEquals(string cdsRef, string cdsAlt, int codonNumRef, int codonNumAlt)
        {
            for (int h = 0, i = 3 * codonNumRef, j = 3 * codonNumAlt; h < 3; i++, j++, h++)
            {
                if (i >= cdsRef.Length || j >= cdsAlt.Length) // Premature end of sequence? (i.e. sequence ends before codon end)
                {
                    return i >= cdsRef.Length && j >= cdsAlt.Length; // We consider them equal only if both sequences reached the end at the same time
                }

                // Same base?
                if (cdsRef[i] != cdsAlt[j])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Get codons from CDS
        /// </summary>
        /// <param name="cds"></param>
        /// <param name="codonNumStart"></param>
        /// <param name="codonNumEnd"></param>
        /// <returns></returns>
        private static string Codons(string cds, int codonNumStart, int codonNumEnd)
        {
            if (codonNumEnd >= 0 && codonNumEnd < codonNumStart) { return ""; }
            int endBase = codonNumEnd < 0 ? cds.Length : 3 * (codonNumEnd + 1);
            endBase = Math.Min(cds.Length, endBase);
            int startBase = Math.Max(0, 3 * codonNumStart);
            int baseLen = endBase - startBase;
            return cds.Substring(startBase, baseLen);
        }

        /// <summary>
        /// Calculate codons by applying the variant and calculating the differences in CDS sequences.
        /// This is a slow method, makes sense only for complex variants
        /// </summary>
        protected void CodonsRefAlt()
        {
            Transcript trNew = Transcript.ApplyVariant(Variant) as Transcript;
            cdsAlt = SequenceExtensions.ConvertToString(trNew.RetrieveCodingSequence());
            cdsRef = SequenceExtensions.ConvertToString(Transcript.RetrieveCodingSequence());
            cdsDiff(); // Calculate differences: CDS
        }

        /// <summary>
        /// How many full / partial exons does the variant affect?
        /// </summary>
        protected void CountAffectedExons()
        {
            exonFull = 0;
            exonPartial = 0;

            foreach (Exon ex in Transcript.Exons)
            {
                if (Variant.Includes(ex)) { exonFull++; }
                else if (Variant.Intersects(ex)) { exonPartial++; }
            }
        }

        protected abstract void EffectTranscript();

        /// <summary>
        /// Variant affect one or more exons
        /// </summary>
        protected abstract void Exons();

        /// <summary>
        /// Variant affect one or more coding exons
        /// </summary>
        protected abstract void ExonsCoding();

        /// <summary>
        /// Variant affect one or more non-coding exons
        /// </summary>
        protected abstract void ExonsNoncoding();

        /// <summary>
        /// Variant affect one intron
        /// </summary>
        protected abstract void Intron();
    }
}