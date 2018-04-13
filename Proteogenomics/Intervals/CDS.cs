using System.Collections.Generic;

namespace Proteogenomics
{
    public class CDS :
        Interval
    {
        public CDS(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
        }

        public override Interval ApplyVariant(Variant variant)
        {
            Interval i = base.ApplyVariant(variant);
            return new CDS(i.Parent as Transcript, i.ChromosomeID, i.Strand, i.OneBasedStart, i.OneBasedEnd, i.Variants);
        }
    }
}