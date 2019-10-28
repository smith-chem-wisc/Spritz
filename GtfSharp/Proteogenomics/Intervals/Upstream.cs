using System.Collections.Generic;

namespace Proteogenomics
{
    public class Upstream
        : Interval
    {
        public Upstream(Transcript parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Upstream(Downstream downstream)
            : base(downstream)
        {
        }
    }
}