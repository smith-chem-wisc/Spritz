namespace Proteogenomics
{
    public class CDS :
        Interval
    {
        public CDS(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public override Interval ApplyVariant(Variant variant)
        {
            IntervalSequence i = base.ApplyVariant(variant) as IntervalSequence;
            return new CDS(i.ChromosomeID, i.Strand, i.OneBasedStart, i.OneBasedEnd);
        }
    }
}