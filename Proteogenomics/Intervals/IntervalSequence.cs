using Bio;
using Bio.Extensions;
using System;
using System.Text;

namespace Proteogenomics
{
    public class IntervalSequence
     : Interval
    {
        #region Public Constructors

        /// <summary>
        /// Sequence using BioDotNet interface
        /// </summary>
        public ISequence Sequence { get; set; }

        public IntervalSequence(Interval parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, ISequence sequence)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd)
        {
            Sequence = sequence;
        }

        public IntervalSequence(Interval interval, ISequence sequence)
            : base(interval)
        {
            Sequence = sequence;
        }

        #endregion Public Constructors

        #region Variant Application Methods

        /// <summary>
        /// Do something with a variant within this interval
        /// </summary>
        /// <param name="variant"></param>
        public override Interval ApplyVariant(Variant variant)
        {
            Variants.Add(variant);

            IntervalSequence newIntervalSeq = new IntervalSequence(base.ApplyVariant(variant), new Sequence(Sequence));
            if (variant.Intersects(this) && Sequence != null && Sequence.Count != 0)
            {
                switch (variant.VarType)
                {
                    case Variant.VariantType.SNV:
                        ApplySnp(variant, newIntervalSeq);
                        break;

                    case Variant.VariantType.INS:
                        ApplyIns(variant, newIntervalSeq);
                        break;

                    case Variant.VariantType.DEL:
                        ApplyDel(variant, newIntervalSeq);
                        break;

                    case Variant.VariantType.DUP:
                        ApplyDup(variant, newIntervalSeq);
                        break;

                    case Variant.VariantType.MNV:
                        ApplyMnp(variant, newIntervalSeq);
                        break;

                    case Variant.VariantType.MIXED:
                        // When applying a MIXED variant, it is decomposed into two
                        // variants (MNP+InDel) and each of them is applied to the
                        // marker. So at this point the variants have been fully
                        // applied and there is no need for further processing.
                        break;

                    default:
                        throw new ArgumentException("Unimplemented method for variant change type " + variant.VarType + "\n\tVariant: " + variant);
                }
            }

            return newIntervalSeq;
        }

        protected void ApplyDel(Variant variant, IntervalSequence markerSeq)
        {
            // Get sequence in positive strand direction
            ISequence seq = IsStrandPlus() ? Sequence : Sequence.GetReverseComplementedSequence();

            // Apply change to sequence
            long idxStart = variant.OneBasedStart - OneBasedStart;
            long idxEnd = idxStart + variant.Length();

            StringBuilder newSeq = new StringBuilder();
            if (idxStart >= 0) { newSeq.Append(SequenceExtensions.ConvertToString(seq, 0, idxStart)); }
            if (idxEnd >= 0 && idxEnd < seq.Count) { newSeq.Append(SequenceExtensions.ConvertToString(seq)); }

            // Update sequence
            seq = new Sequence(seq.Alphabet, newSeq.ToString());
            markerSeq.Sequence = IsStrandPlus() ? seq : seq.GetReverseComplementedSequence();
        }

        protected void ApplyDup(Variant variant, IntervalSequence markerSeq)
        {
            // Get sequence in positive strand direction
            ISequence seq = IsStrandPlus() ? Sequence : Sequence.GetReverseComplementedSequence();

            // Apply duplication to sequence
            String dupSeq = SequenceExtensions.ConvertToString(GetSequence(Intersect(variant)));
            long idx = variant.OneBasedStart - OneBasedStart - 1;
            if (idx >= 0)
            {
                seq = new Sequence(seq.Alphabet, SequenceExtensions.ConvertToString(seq, 0, idx + 1) + dupSeq + SequenceExtensions.ConvertToString(seq, idx + 1));
            }
            else
            {
                seq = new Sequence(seq.Alphabet, dupSeq + SequenceExtensions.ConvertToString(seq));
            }

            // Update sequence
            markerSeq.Sequence = IsStrandPlus() ? seq : seq.GetReverseComplementedSequence();
        }

        protected void ApplyIns(Variant variant, IntervalSequence markerSeq)
        {
            // Get sequence in positive strand direction
            ISequence seq = IsStrandPlus() ? Sequence : Sequence.GetReverseComplementedSequence();

            // Apply change to sequence
            String netChange = variant.NetChange(this);
            long idx = variant.OneBasedStart - OneBasedStart - 1;
            if (idx >= 0)
            {
                seq = new Sequence(seq.Alphabet, SequenceExtensions.ConvertToString(seq, 0, idx + 1) + netChange + SequenceExtensions.ConvertToString(seq, idx + 1));
            }
            else
            {
                seq = new Sequence(seq.Alphabet, netChange + SequenceExtensions.ConvertToString(seq));
            }

            // Update sequence
            markerSeq.Sequence = IsStrandPlus() ? seq : seq.GetReverseComplementedSequence();
        }

        protected void ApplyMnp(Variant variant, IntervalSequence markerSeq)
        {
            // Calculate indexes
            long idxStart = variant.OneBasedStart - OneBasedStart;
            long idxAlt = 0;

            // Variant starts before this marker (e.g. motif with sequence)
            if (idxStart < 0)
            {
                idxAlt = -idxStart; // Remove first 'idxStart' bases from ALT sequence
                idxStart = 0;
            }

            long changeSize = variant.IntersectSize(this);
            long idxEnd = idxStart + changeSize;

            // Apply variant to sequence
            ISequence seq = IsStrandPlus() ? Sequence : Sequence.GetReverseComplementedSequence();
            StringBuilder seqsb = new StringBuilder();
            seqsb.Append(SequenceExtensions.ConvertToString(seq, 0, idxStart));
            String seqAlt = variant.SecondAlleleString.Substring((int)idxAlt, (int)(idxAlt + changeSize));
            seqsb.Append(seqAlt);
            seqsb.Append(SequenceExtensions.ConvertToString(seq, idxEnd));

            // Update sequence
            seq = new Sequence(seq.Alphabet, seqsb.ToString());
            markerSeq.Sequence = IsStrandPlus() ? seq : seq.GetReverseComplementedSequence();
        }

        protected void ApplySnp(Variant variant, IntervalSequence markerSeq)
        {
            // Get sequence in positive strand direction
            ISequence seq = IsStrandPlus() ? Sequence : Sequence.GetReverseComplementedSequence();

            // Apply change to sequence
            long idx = variant.OneBasedStart - OneBasedStart;
            seq = new Sequence(seq.Alphabet, SequenceExtensions.ConvertToString(seq, 0, idx) + variant.SecondAlleleString + SequenceExtensions.ConvertToString(seq, idx + 1));

            // Update sequence
            markerSeq.Sequence = IsStrandPlus() ? seq : seq.GetReverseComplementedSequence();
        }

        #endregion Variant Application Methods

        #region Other Methods

        /// <summary>
        /// ase in this marker at position 'index' (relative to marker start)
        /// </summary>
        /// <param name="index"></param>
        /// <param name="len"></param>
        /// <returns></returns>
        public ISequence basesAt(long index, long len)
        {
            if (IsStrandMinus())
            {
                long idx = Sequence.Count - index - len;
                return Sequence.GetSubSequence(idx, len).GetReverseComplementedSequence(); // Minus strand => Sequence has been reversed and WC-complemented
            }

            return Sequence.GetSubSequence(index, len);
        }

        /// <summary>
        /// Base at position 'pos' (genomic coordinates)
        /// </summary>
        /// <param name="pos"></param>
        /// <param name="len"></param>
        /// <returns></returns>
        public ISequence basesAtPos(int pos, int len)
        {
            long index = pos - OneBasedEnd;
            if (index < 0) { return new Sequence(Alphabets.DNA, ""); }
            return basesAt(index, len);
        }

        public ISequence GetSequence(Interval interval)
        {
            if (!Intersects(interval)) { return null; }
            return Sequence.GetSubSequence(interval.OneBasedStart, interval.Length());
        }

        #endregion Other Methods
    }
}