using System.Collections.Generic;

namespace Proteogenomics
{
    public class UTR3Prime
        : UTR
    {
        public UTR3Prime(Exon parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public UTR3Prime(UTR3Prime utr)
            : base(utr)
        {
        }

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "three_prime_UTR";

        public override bool is3Prime()
        {
            return true;
        }

        public override bool is5Prime()
        {
            return false;
        }
    }
}