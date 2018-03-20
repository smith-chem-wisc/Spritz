using Bio;
using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class Exon :
        IntervalSequence
    {
        public Exon(ISequence Sequence, long oneBasedStart, long oneBasedEnd, string chromID, string strand)
            : base(chromID, strand, oneBasedStart, oneBasedEnd, Sequence)
        {
        }

        public override Interval ApplyVariant(Variant variant)
        {
            IntervalSequence i = base.ApplyVariant(variant) as IntervalSequence;
            return new Exon(i.Sequence, i.OneBasedStart, i.OneBasedEnd, i.ChromosomeID, i.Strand);
        }
    }
}