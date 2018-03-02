using Bio;

namespace Proteogenomics
{
    public class Exon :
        IntervalSequence
    {
        public Exon(ISequence Sequence, long oneBasedStart, long oneBasedEnd, string chromID, string strand)
            : base(chromID, strand, oneBasedStart, oneBasedEnd, Sequence)
        {
        }
    }
}