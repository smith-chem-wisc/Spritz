namespace Proteogenomics
{
    public class Upstream
        : Interval
    {
        public Upstream(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Upstream(Downstream downstream)
            : base(downstream)
        {
        }
    }
}