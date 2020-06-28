using System.Collections.Generic;

namespace Proteogenomics
{
    public class Intron
        : Interval
    {
        public Intron(Transcript parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Intron(Intron intron)
            : base(intron)
        {
        }
    }
}