namespace Proteogenomics
{
    public class Downstream
        : Interval
    {
        public Downstream(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Downstream(Downstream downstream)
            : base(downstream)
        {
        }
    }
}