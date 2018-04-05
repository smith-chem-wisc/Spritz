using Bio;
using Bio.Extensions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class UTR5Prime
        : UTR
    {
        private List<UTR5Prime> UTRs { get; set; }

        public UTR5Prime(Exon parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants) :
            base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
        }

        public UTR5Prime(UTR5Prime utr) : base(utr)
        {
        }

        public override bool is3Prime()
        {
            return false;
        }

        public override bool is5Prime()
        {
            return true;
        }

        private long utrDistance(Variant variant, Transcript tr)
        {
            long cdsStart = tr.cdsStart;
            if (cdsStart < 0)
            {
                return -1;
            }
            if (IsStrandPlus())
            {
                return cdsStart - variant.OneBasedEnd;
            }
            return variant.OneBasedStart - cdsStart;
        }

        public List<UTR5Prime> get5primeUtrs()
        {
            if (UTRs == null)
            {
                Transcript tr = (Transcript)FindParent(typeof(Transcript));

                // Get UTRs and sort them
                UTRs = tr.UTRs.OfType<UTR5Prime>().ToList();
                if (IsStrandPlus())
                {
                    UTRs = UTRs.OrderBy(u => u.OneBasedStart).ToList(); // Sort by start position
                }
                else
                {
                    UTRs = UTRs.OrderBy(u => u.OneBasedEnd).ToList(); // Sort by end position (reversed)
                }
            }

            return UTRs;
        }

        public string getSequence()
        {
            // Create UTR sequence
            StringBuilder sb = new StringBuilder();
            foreach (UTR5Prime utr in get5primeUtrs())
            {
                Exon ex = (Exon)utr.Parent;
                ISequence utrSeq = ex.Sequence;
                if (utr.Length() < utrSeq.Count) { utrSeq = utrSeq.GetSubSequence(0, utr.Length()); } // UTR5' may stop before end of exon
                sb.Append(SequenceExtensions.ConvertToString(utrSeq));
            }
            return sb.ToString();
        }

        /// <summary>
        /// Is a new start codon produced? @return New start codon (or empty string if there is no new start codon)
        /// </summary>
        /// <param name="chars"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        private static string startGained(char[] chars, long pos)
        {
            // Analyze all frames
            for (long i = Math.Max(0, pos - 2); (i <= pos) && ((i + 2) < chars.Length); i++)
            {
                string codon = "" + chars[i] + chars[i + 1] + chars[i + 2];
                if (CodonsStandard.START_CODONS.Contains(codon.ToUpper(CultureInfo.InvariantCulture)))
                {
                    return codon.ToUpper(CultureInfo.InvariantCulture); // This frame has a start codon?
                }
            }
            return "";
        }

        /// <summary>
        /// Did we gain a start codon in this 5'UTR interval?
        /// </summary>
        /// <param name="seqChange"></param>
        /// <returns>A new start codon (if gained)</returns>
        private string startGained(Variant seqChange)
        {
            if (!seqChange.isSnv()) { return ""; } // Only SNPs supported.

            // Calculate SNP position relative to UTRs
            long pos = seqChange.DistanceBases(get5primeUtrs().OfType<Interval>().ToList(), IsStrandMinus());

            // Change base at SNP position
            string sequence = getSequence();
            char[] chars = sequence.ToCharArray();
            char snpBase = seqChange.NetChange(this)[0];
            if (IsStrandMinus() && Alphabets.DNA.TryGetComplementSymbol((byte)snpBase, out byte complement))
            {
                snpBase = (char)complement;
            }
            chars[pos] = snpBase;

            // Do we gain a new start codon?
            return startGained(chars, pos);
        }

        public override bool CreateVariantEffect(Variant variant, VariantEffects variantEffects)
        {
            // Has the whole UTR been deleted?
            if (variant.Includes(this) && (variant.VarType == Variant.VariantType.DEL))
            {
                variantEffects.add(variant, this, EffectType.UTR_5_DELETED, ""); // A UTR was removed entirely
                return true;
            }

            // Add distance
            Transcript tr = (Transcript)FindParent(typeof(Transcript));
            long distance = utrDistance(variant, tr);
            VariantEffect variantEffect = new VariantEffect(variant);
            variantEffect.set(this, IntervalType, VariantEffect.EffectDictionary[IntervalType], distance >= 0 ? distance + " bases from TSS" : "");
            variantEffect.setDistance(distance);
            variantEffects.add(variantEffect);

            // Start gained?
            string gained = startGained(variant);
            if (gained != "")
            {
                variantEffects.add(variant, this, EffectType.START_GAINED, gained);
            }

            return true;
        }
    }
}