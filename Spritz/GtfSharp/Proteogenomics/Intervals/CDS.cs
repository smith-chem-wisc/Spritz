using System.Collections.Generic;

namespace Proteogenomics
{
    public class CDS :
        Interval
    {
        public CDS(Transcript parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd, int startFrame)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
            StartFrame = startFrame;
        }

        public CDS(CDS cds)
            : base(cds)
        {
        }

        /// <summary>
        /// Start frame for this coding domain sequence { 0, 1, 2 } and -1 is undefined, although that isn't implemented in this program
        /// </summary>
        public int StartFrame { get; private set; }

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "CDS";

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
                utr5 = new UTR5Prime(parent, parent.Source, parent.ChromosomeID, parent.Strand, OneBasedStart, end);
            }
            else
            {
                long start = OneBasedEnd - (StartFrame - 1);
                utr5 = new UTR5Prime(parent, parent.Source, parent.ChromosomeID, parent.Strand, start, OneBasedEnd);
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
                utr3 = new UTR3Prime(parent, parent.ChromosomeID, parent.Source, parent.Strand, start, OneBasedEnd);
            }
            else
            {
                long end = OneBasedStart + (endFrame - 1);
                utr3 = new UTR3Prime(parent, parent.ChromosomeID, parent.Source, parent.Strand, OneBasedStart, end);
            }

            // correct start or end coordinates
            if (IsStrandPlus()) { OneBasedEnd -= endFrame; }
            else { OneBasedStart += endFrame; }

            return utr3;
        }
    }
}