using System.Collections.Generic;

namespace Proteogenomics
{
    public class CDS :
        Interval
    {
        public CDS(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants, int startFrame)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
            StartFrame = startFrame;
        }

        public CDS(CDS cds)
            : base(cds)
        {
        }

        public int StartFrame { get; private set; }

        /// <summary>
        /// Corrrects coordinates of end of coding sequence based on coding sequence length
        /// </summary>
        /// <param name="codingSequenceLength"></param>
        /// <returns></returns>
        public UTR5Prime StartFrameCorrection()
        {
            if (StartFrame <= 0) { return null; } // nothing to do

            // First exon is not zero? => Create a UTR5 prime to compensate
            Exon parent = (FindParent(typeof(Transcript)) as Transcript).GetLastCodingExon();
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
            if (IsStrandPlus()) { OneBasedStart += StartFrame; }
            else { OneBasedEnd -= StartFrame; }

            StartFrame = 0;

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
            return new CDS(i.Parent as Transcript, i.ChromosomeID, i.Strand, i.OneBasedStart, i.OneBasedEnd, i.Variants, this.StartFrame);
        }
    }
}