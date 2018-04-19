using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class CDS :
        Interval
    {
        private static HashSet<int> possibleFrames = new HashSet<int>{ -1, 0, 1, 2 };

        private int _StartFrame;

        public CDS(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, int frame, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
            StartFrame = frame;
        }

        /// <summary>
        /// Frame can be {-1, 0, 1, 2}, where -1 means unknown
        /// </summary>
        public int StartFrame
        {
            get
            {
                return _StartFrame;
            }

            private set
            {
                if (!possibleFrames.Contains(value)) { throw new ArgumentOutOfRangeException("frame"); }
                else { _StartFrame = value; }
            }
        }

        /// <summary>
        /// Correct coordinates according to frame differences
        /// </summary>
        /// <param name="correction"></param>
        /// <returns></returns>
        public UTR5Prime StartFrameCorrection(int correction)
        {
            if (correction <= 0) { return null; } // nothing to do
            if (Length() <= correction) { return null; } // too short, cannot correct frame

            // First exon is not zero? => Create a UTR5 prime to compensate
            Exon parent = (FindParent(typeof(Transcript)) as Transcript).GetFirstCodingExon();

            UTR5Prime utr5 = null;
            if (IsStrandPlus())
            {
                long end = OneBasedStart + (StartFrame - 1);
                utr5 = new UTR5Prime(parent, parent.ChromosomeID, parent.Strand, OneBasedStart, end, null);
            }
            else
            {
                long start = OneBasedEnd - (StartFrame - 1);
                utr5 = new UTR5Prime(parent, parent.ChromosomeID, parent.Strand, start, OneBasedEnd, null);
            }

            // correct start or end coordinates
            if (IsStrandPlus()) { OneBasedStart += correction; }
            else { OneBasedEnd -= correction; }

            // correct frame
            StartFrame = (StartFrame - correction) % 3;
            while (StartFrame < 0)
            {
                StartFrame += 3;
            }

            return utr5;
        }

        /// <summary>
        /// Corrrects coordinates of end of coding sequence based on coding sequence length
        /// </summary>
        /// <param name="codingSequenceLength"></param>
        /// <returns></returns>
        public UTR3Prime EndFrameCorrection(long codingSequenceLength)
        {
            long endFrame = codingSequenceLength % 3;
            if (endFrame <= 0) { return null; } // nothing to do

            // First exon is not zero? => Create a UTR5 prime to compensate
            Exon parent = (FindParent(typeof(Transcript)) as Transcript).GetLastCodingExon();
            UTR3Prime utr3 = null;
            if (IsStrandPlus())
            {
                long start = OneBasedEnd - (endFrame - 1);
                utr3 = new UTR3Prime(parent, parent.ChromosomeID, parent.Strand, start, OneBasedEnd, null);
            }
            else
            {
                long end = OneBasedStart + (endFrame - 1);
                utr3 = new UTR3Prime(parent, parent.ChromosomeID, parent.Strand, OneBasedStart, end, null);
            }

            // correct start or end coordinates
            if (IsStrandPlus()) { OneBasedEnd -= endFrame; }
            else { OneBasedStart += endFrame; }

            return utr3;
        }

        public override Interval ApplyVariant(Variant variant)
        {
            Interval i = base.ApplyVariant(variant);
            return new CDS(i.Parent as Transcript, i.ChromosomeID, i.Strand, i.OneBasedStart, i.OneBasedEnd, this.StartFrame, i.Variants);
        }
    }
}