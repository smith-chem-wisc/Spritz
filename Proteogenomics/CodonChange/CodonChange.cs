using Bio;
using Bio.Extensions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace Proteogenomics
{
    public class CodonChange
    {
        public static readonly HashSet<string> StartCodons = new HashSet<string> { "ATG" };
        public static readonly bool ShowCodonChange = true; // This is disabled in some specific test cases
        public static readonly int CODON_SIZE = 3; // I'll be extremely surprised if you ever need to change this parameter...

        public int CodonStartNumber { get; set; } = -1;
        public int CodonStartIndex { get; set; } = -1;
        protected bool ReturnNow { get; set; } // Can we return immediately after calculating the first 'codonChangeSingle()'?
        protected bool RequireNetCdsChange { get; set; }
        protected Variant Variant { get; set; }
        protected Transcript Transcript { get; set; }
        protected Exon Exon { get; set; }
        protected VariantEffects VariantEffects { get; set; }
        protected string CodonsReference { get; set; } = ""; // REF codons (without variant)
        protected string CodonsAlternate { get; set; } = ""; // ALT codons (after variant is applied)
        protected string NetCodingSequenceChange { get; set; } = "";

        /// <summary>
        /// Create a specific codon change for a variant
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="transcript"></param>
        /// <param name="variantEffects"></param>
        /// <returns></returns>
        public static CodonChange Factory(Variant variant, Transcript transcript, VariantEffects variantEffects)
        {
            switch (variant.VarType)
            {
                case Variant.VariantType.SNV:
                    return new CodonChangeSnv(variant, transcript, variantEffects);

                case Variant.VariantType.INS:
                    return new CodonChangeIns(variant, transcript, variantEffects);

                case Variant.VariantType.DEL:
                    return new CodonChangeDel(variant, transcript, variantEffects);

                case Variant.VariantType.MNV:
                    return new CodonChangeMnv(variant, transcript, variantEffects);

                //case Variant.VariantType.MIXED:
                //    return new CodonChangeMixed(variant, transcript, variantEffects);

                case Variant.VariantType.DUP:
                    return new CodonChangeDup(variant, transcript, variantEffects);

                case Variant.VariantType.INV:
                    return new CodonChangeInv(variant, transcript, variantEffects);

                case Variant.VariantType.INTERVAL:
                    return new CodonChangeInterval(variant, transcript, variantEffects);

                default:
                    throw new ArgumentException("Unimplemented factory for variant type '" + variant.VarType + "', variant: " + variant);
            }
        }

        protected CodonChange(Variant variant, Transcript transcript, VariantEffects variantEffects)
        {
            Transcript = transcript;
            VariantEffects = variantEffects;
            Variant = variant;
        }

        /// <summary>
        /// Calculate additional effect due to codon changes, e.g. A frame-shift that also affects a stop codon
        /// </summary>
        /// <param name="codonsOld"></param>
        /// <param name="codonsNew"></param>
        /// <param name="codonNum"></param>
        /// <param name="codonIndex"></param>
        /// <param name="aaOld"></param>
        /// <param name="aaNew"></param>
        /// <returns></returns>
        protected EffectType AdditionalEffect(string codonsOld, string codonsNew, int codonNum, int codonIndex, string aaOld, string aaNew)
        {
            EffectType newEffectType = EffectType.NONE;

            bool hasOldAa = CodonExtensions.TryTranslateCodon(Transcript.Gene.Chromosome.Mitochondrial, codonsOld, out byte oldAA);
            bool hasNewAa = CodonExtensions.TryTranslateCodon(Transcript.Gene.Chromosome.Mitochondrial, codonsNew, out byte newAA);
            bool isStartNew = Transcript.Gene.Chromosome.Mitochondrial ?
                CodonsVertebrateMitochondrial.START_CODONS.Contains(codonsNew.ToUpper(CultureInfo.InvariantCulture)) :
                CodonsStandard.START_CODONS.Contains(codonsNew.ToUpper(CultureInfo.InvariantCulture));
            bool isStartOld = Transcript.Gene.Chromosome.Mitochondrial ?
                CodonsVertebrateMitochondrial.START_CODONS.Contains(codonsOld.ToUpper(CultureInfo.InvariantCulture)) :
                CodonsStandard.START_CODONS.Contains(codonsOld.ToUpper(CultureInfo.InvariantCulture));
            bool isStopOld = hasOldAa && oldAA == Alphabets.Protein.Ter;
            bool isStopNew = hasNewAa && newAA == Alphabets.Protein.Ter;

            if (Variant.isSnv() || Variant.isMnv())
            {
                // SNM and MNP effects
                if (aaOld.Equals(aaNew))
                {
                    // Same AA: Synonymous coding
                    if (codonNum == 0 && isStartOld)
                    {
                        // It is in the first codon (which also is a start codon)
                        if (isStartNew) { newEffectType = EffectType.SYNONYMOUS_START; }// The new codon is also a start codon => SYNONYMOUS_START
                        else { newEffectType = EffectType.START_LOST; } // The AA is the same, but the codon is not a start codon => start lost
                    }
                    else if (isStopOld)
                    {
                        // Stop codon
                        if (isStopNew) { newEffectType = EffectType.SYNONYMOUS_STOP; }// New codon is also a stop => SYNONYMOUS_STOP
                        else { newEffectType = EffectType.STOP_LOST; }// New codon is not a stop, the we've lost a stop
                    }
                    else
                    {
                        newEffectType = EffectType.SYNONYMOUS_CODING;  // All other cases are just SYNONYMOUS_CODING
                    }
                }
                else
                {
                    // Different AA: Non-synonymous coding
                    if ((codonNum == 0) && isStartOld)
                    {
                        // It is in the first codon (which also is a start codon)
                        if (isStartNew) { newEffectType = EffectType.NON_SYNONYMOUS_START; }// Non-synonymous mutation on first codon => start lost
                        else { newEffectType = EffectType.START_LOST; }// Non-synonymous mutation on first codon => start lost
                    }
                    else if (isStopOld)
                    {
                        // Stop codon
                        if (isStopNew) { newEffectType = EffectType.NON_SYNONYMOUS_STOP; } // Notice: This should never happen for SNPs! (for some reason I removed this comment at some point and that create some confusion): http://www.biostars.org/post/show/51352/in-snpeff-impact-what-is-difference-between-stop_gained-and-non-synonymous_stop/
                        else { newEffectType = EffectType.STOP_LOST; }
                    }
                    else if (isStopNew)
                    {
                        newEffectType = EffectType.STOP_GAINED;
                    }
                    else
                    {
                        newEffectType = EffectType.NON_SYNONYMOUS_CODING; // All other cases are just NON_SYN
                    }
                }
            }
            else
            {
                // Add a new effect in some cases
                if ((codonNum == 0) && isStartOld && !isStartNew)
                {
                    newEffectType = EffectType.START_LOST;
                }
                else if (isStopOld && !isStopNew)
                {
                    newEffectType = EffectType.STOP_LOST;
                }
                else if (!isStopOld && isStopNew)
                {
                    newEffectType = EffectType.STOP_GAINED;
                }
            }

            return newEffectType;
        }

        /// <summary>
        /// Calculate base number in a cds where 'pos' is
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        protected long CdsBaseNumber(long pos)
        {
            long cdsbn = Transcript.BaseNumberCds(pos, true);

            // Does not intersect the transcript?
            if (cdsbn < 0)
            {
                // 'pos' before transcript start
                if (pos <= Transcript.cdsStart)
                {
                    if (Transcript.IsStrandPlus()) { return 0; }
                    return Transcript.RetrieveCodingSequence().Count;
                }

                // 'pos' is after CDS end
                if (Transcript.IsStrandPlus()) { return Transcript.RetrieveCodingSequence().Count; }
                return 0;
            }

            return cdsbn;
        }

        /// <summary>
        /// Calculate a list of codon changes
        /// </summary>
        public virtual void ChangeCodon()
        {
            if (!Transcript.Intersects(Variant))
            {
                return;
            }

            // Get coding start (after 5 prime UTR)
            long cdsStart = Transcript.cdsStart;

            // We may have to calculate 'netCdsChange', which is the effect on the CDS
            NetCodingSequenceChange = NetCdsChange();
            if (RequireNetCdsChange && NetCodingSequenceChange == "")
            { // This can happen on mixed changes where the 'InDel' part lies outside the transcript's exons
                CodonsReference = "";
                CodonsAlternate = "";
                return;
            }

            //---
            // Concatenate all exons
            //---
            int firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
            List<Exon> exons = Transcript.ExonsSortedStrand;
            foreach (Exon exon in exons)
            {
                Exon = exon;
                if (exon.Intersects(Variant))
                {
                    long cdsBaseInExon; // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)

                    if (Transcript.IsStrandPlus())
                    {
                        long firstvariantBaseInExon = Math.Max(Variant.OneBasedStart, Math.Max(exon.OneBasedStart, cdsStart));
                        cdsBaseInExon = firstvariantBaseInExon - Math.Max(exon.OneBasedStart, cdsStart);
                    }
                    else
                    {
                        long lastvariantBaseInExon = Math.Min(Variant.OneBasedEnd, Math.Min(exon.OneBasedEnd, cdsStart));
                        cdsBaseInExon = Math.Min(exon.OneBasedEnd, cdsStart) - lastvariantBaseInExon;
                    }

                    if (cdsBaseInExon < 0) { cdsBaseInExon = 0; }

                    // Get codon number and index within codon (where seqChage is pointing)
                    if (CodonStartNumber < 0)
                    {
                        CodonStartNumber = (int)(firstCdsBaseInExon + cdsBaseInExon) / CODON_SIZE;
                        CodonStartIndex = (int)(firstCdsBaseInExon + cdsBaseInExon) % CODON_SIZE;
                    }

                    // Use appropriate method to calculate codon change
                    //bool hasChanged = false; // Was there any change?
                    //hasChanged = ChangeCodon(exon);

                    // Any change? => Add change to list
                    //if (hasChanged && !VariantEffects.hasMarker()) VariantEffects.setMarker(exon); // It is affecting this exon, so we set the marker

                    // Can we finish after effect of first exon is added?
                    if (ReturnNow) { return; }
                }

                firstCdsBaseInExon += Transcript.IsStrandPlus() ?
                    (int)Math.Max(0, exon.OneBasedEnd - Math.Max(exon.OneBasedStart, cdsStart) + 1) :
                    (int)Math.Max(0, Math.Min(cdsStart, exon.OneBasedEnd) - exon.OneBasedStart + 1);
            }
        }

        /// <summary>
        /// Calculate the effect on an exon
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        protected virtual bool ChangeCodon(Exon exon)
        {
            throw new InvalidOperationException("Unimplemented method codonChangeSingle() for\n\t\tVariant type : " + Variant.GetType().Name + "\n\t\tClass        : " + GetType().Name + "\n\t\tVariant      : " + Variant);
        }

        /// <summary>
        /// Calculate new codons
        /// </summary>
        /// <returns></returns>
        protected virtual string CodonsAlt()
        {
            throw new InvalidOperationException("Unimplemented method for this thype of CodonChange: " + this.GetType().Name);
        }

        /// <summary>
        /// Calculate 'reference' codons
        /// </summary>
        /// <returns></returns>
        protected virtual string CodonsRef()
        {
            return CodonsRef(1);
        }

        /// <summary>
        /// Calculate 'reference' codons
        /// </summary>
        /// <param name="numCodons"></param>
        /// <returns></returns>
        protected string CodonsRef(int numCodons)
        {
            ISequence cds = Transcript.RetrieveCodingSequence();
            string codon = "";

            int start = CodonStartNumber * CODON_SIZE;
            int end = start + numCodons * CODON_SIZE;

            int len = (int)cds.Count;
            if (start >= len) { start = len; }
            if (end >= len) { end = len; }

            // Capitalize
            codon = SequenceExtensions.ConvertToString(cds.GetSubSequence(start, end));

            // Codon not multiple of three? Add missing bases as 'N'
            if (codon.Length % 3 == 1) { codon += "NN"; }
            else if (codon.Length % 3 == 2) { codon += "N"; }

            return codon;
        }

        protected VariantEffect Effect(Interval marker, EffectType effectType, bool allowReplace)
        {
            return Effect(marker, effectType, VariantEffect.EffectDictionary[effectType], CodonsReference, CodonsAlternate, CodonStartNumber, CodonStartIndex, allowReplace);
        }

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        /// <param name="effectImpact"></param>
        /// <param name="codonsOld"></param>
        /// <param name="codonsNew"></param>
        /// <param name="codonNum"></param>
        /// <param name="codonIndex"></param>
        /// <param name="allowReplace"></param>
        /// <returns></returns>
        private VariantEffect Effect(Interval marker, EffectType effectType, EffectImpact effectImpact, string codonsOld, string codonsNew, int codonNum, int codonIndex, bool allowReplace)
        {
            // Create and add variant affect
            long cDnaPos = Transcript.BaseNumber2MRnaPos(Variant.OneBasedStart);
            VariantEffect varEff = new VariantEffect(Variant, marker, effectType, effectImpact, codonsOld, codonsNew, codonNum, codonIndex, cDnaPos);
            VariantEffects.add(varEff);

            // Are there any additional effects? Sometimes a new effect arises from setting codons (e.g. FRAME_SHIFT disrupts a STOP codon)
            EffectType addEffType = AdditionalEffect(codonsOld, codonsNew, codonNum, codonIndex, varEff.aaRef, varEff.aaAlt);
            if (addEffType != EffectType.NONE && addEffType != effectType)
            {
                if (allowReplace && addEffType.CompareTo(effectType) < 0)
                {
                    // Replace main effect (using default impact)
                    varEff.setEffect(addEffType);
                }
                else
                {
                    // Add effect to list (using default impact)
                    varEff.addEffect(addEffType);
                }
            }

            return varEff;
        }

        protected VariantEffect EffectNoCodon(Interval marker, EffectType effectType)
        {
            return Effect(marker, effectType, VariantEffect.EffectDictionary[effectType], "", "", -1, -1, false);
        }

        protected VariantEffect EffectNoCodon(Interval marker, EffectType effectType, EffectImpact effectImpact)
        {
            return Effect(marker, effectType, effectImpact, "", "", -1, -1, false);
        }

        /// <summary>
        /// Does the variant intersect any exons?
        /// </summary>
        /// <returns></returns>
        protected bool IntersectsExons()
        {
            foreach (Exon ex in Transcript.Exons)
            {
                if (Variant.Intersects(ex))
                {
                    return true;
                }
            }
            return false;
        }

        /// <summary>
        /// We may have to calculate 'netCdsChange', which is the effect on the CDS.
        /// Note: A deletion or a MNP might affect several exons
        /// </summary>
        /// <returns></returns>
        protected virtual string NetCdsChange()
        {
            if (!RequireNetCdsChange) { return ""; }

            if (Variant.Length() > 1)
            {
                StringBuilder sb = new StringBuilder();
                foreach (Exon exon in Transcript.ExonsSortedStrand)
                {
                    sb.Append(Variant.NetChange(exon));
                }
                return sb.ToString();
            }

            return Variant.NetChange(Transcript.IsStrandMinus());
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append("Transcript : " + Transcript.ID + "\n");
            sb.Append("Variant    : " + Variant + "\n");
            sb.Append("Codons     : " + CodonsReference + "/" + CodonsAlternate + "\tnum: " + CodonStartNumber + "\tidx: " + CodonStartIndex + "\n");
            sb.Append("Effects    :\n");
            foreach (VariantEffect veff in VariantEffects.Effects)
            {
                sb.Append("\t" + veff.getEffectTypeString(false) + "\t" + veff.codonsRef + "/" + veff.codonsAlt + "\t" + veff.aaRef + "/" + veff.aaAlt + "\n");
            }

            return sb.ToString();
        }
    }
}