using Bio;
using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class Exon :
        IntervalSequence
    {
        /// <summary>
        /// Construct an exon
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="Sequence"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chromID"></param>
        /// <param name="strand"></param>
        public Exon(Transcript parent, ISequence Sequence, long oneBasedStart, long oneBasedEnd, string chromID, string strand, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, Sequence, variants)
        {
        }

        /// <summary>
        /// Copy an exon
        /// </summary>
        /// <param name="x"></param>
        public Exon(Exon x)
            : this(x.Parent as Transcript, x.Sequence, x.OneBasedStart, x.OneBasedEnd, x.ChromosomeID, x.Strand, x.Variants)
        {
        }

        /// <summary>
        /// Apply a variant to this exon interval and sequence
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        public override Interval ApplyVariant(Variant variant)
        {
            IntervalSequence i = base.ApplyVariant(variant) as IntervalSequence;
            return new Exon(i.Parent as Transcript, i.Sequence, i.OneBasedStart, i.OneBasedEnd, i.ChromosomeID, i.Strand, i.Variants);
        }

        /// <summary>
        /// Check that the base in the exon corresponds with the one in the SNP
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        public ErrorWarningType SanityCheck(Variant variant)
        {
            if (!Intersects(variant)) { return ErrorWarningType.NONE; }

            // Only makes sense for SNPs and MNPs
            if ((variant.VarType != Variant.VariantType.SNV) && (variant.VarType != Variant.VariantType.MNV)) return ErrorWarningType.NONE;

            // Calculate reference sequence coordinates
            long mstart = Math.Max(variant.OneBasedStart, OneBasedStart);
            long idxStart = mstart - OneBasedStart;

            if (Sequence.Count <= 0) { return ErrorWarningType.WARNING_SEQUENCE_NOT_AVAILABLE; }
            if (idxStart >= Sequence.Count) { return ErrorWarningType.ERROR_OUT_OF_EXON; }

            long mend = Math.Min(variant.OneBasedEnd, OneBasedEnd);
            long len = mend - mstart + 1;

            ISequence realReference = basesAt(idxStart, len);

            // Get variant's reference coordinates
            long varRefStart = mstart - variant.OneBasedStart;
            if (varRefStart < 0) { return ErrorWarningType.ERROR_OUT_OF_EXON; }

            long varRefEnd = mend - variant.OneBasedStart;

            // Variant's reference sequence
            string refStr;
            //if (variant.isNonRef()) // used for some advanced stuff in SnpEff, comparing cancer vs somatic tissue
            //{
            //    refStr = ((VariantNonRef)variant).getVariantRef().getAlt();
            //}
            //else
            //{
            refStr = variant.ReferenceAlleleString;
            //}

            if (varRefEnd >= refStr.Length) { return ErrorWarningType.ERROR_OUT_OF_EXON; }

            int varRefLen = (int)(varRefEnd - varRefStart + 1);
            string variantReference = refStr.Substring((int)varRefStart, varRefLen);

            // Reference sequence different than expected?
            if (!realReference.Equals(variantReference))
            { //
                return ErrorWarningType.WARNING_REF_DOES_NOT_MATCH_GENOME;
            }

            // OK
            return ErrorWarningType.NONE;
        }
    }
}