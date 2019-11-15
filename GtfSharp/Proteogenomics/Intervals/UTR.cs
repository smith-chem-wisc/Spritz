using System.Collections.Generic;

namespace Proteogenomics
{
    public abstract class UTR :
        Interval
    {
        protected UTR(Exon parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        protected UTR(UTR utr)
            : base(utr)
        {
        }

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "UTR";

        public abstract bool is3Prime();

        public abstract bool is5Prime();
    }
}