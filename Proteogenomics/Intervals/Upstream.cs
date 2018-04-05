using System.Collections.Generic;

namespace Proteogenomics
{
    public class Upstream
        : Interval
    {
        public Upstream(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
        }

        public Upstream(Downstream downstream)
            : base(downstream)
        {
        }
    }
}