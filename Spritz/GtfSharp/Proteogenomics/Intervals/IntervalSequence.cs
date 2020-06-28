using Bio;
using Bio.Extensions;
using System;
using System.Collections.Generic;
using System.Text;

namespace Proteogenomics
{
    public class IntervalSequence
     : Interval
    {
        public IntervalSequence(Interval parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd, ISequence sequence)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
            Sequence = sequence;
        }

        public IntervalSequence(Interval interval, ISequence sequence)
            : base(interval)
        {
            Sequence = sequence;
        }

        /// <summary>
        /// Sequence using BioDotNet interface
        /// </summary>
        public ISequence Sequence { get; set; }

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