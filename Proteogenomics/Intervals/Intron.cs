namespace Proteogenomics
{
    public class Intron
        : Interval
    {
        public Intron(Transcript parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Intron(Intron intron)
            : base(intron)
        {
        }
    }
}