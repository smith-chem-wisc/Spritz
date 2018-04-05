using System.Collections.Generic;

namespace Proteogenomics
{
    public class Intron
        : Interval
    {
        public Intron(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
        }

        public Intron(Intron intron)
            : base(intron)
        {
        }
    }
}