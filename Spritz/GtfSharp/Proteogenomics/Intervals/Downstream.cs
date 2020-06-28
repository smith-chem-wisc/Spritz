using System.Collections.Generic;

namespace Proteogenomics
{
    public class Downstream
        : Interval
    {
        public Downstream(Transcript parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Downstream(Downstream downstream)
            : base(downstream)
        {
        }
    }
}