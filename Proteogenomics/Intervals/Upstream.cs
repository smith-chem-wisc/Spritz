namespace Proteogenomics
{
    public class Upstream
        : Interval
    {
        #region Constructors

        public Upstream(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Upstream(Downstream downstream)
            : base(downstream)
        {
        }

        #endregion Constructors
    }
}