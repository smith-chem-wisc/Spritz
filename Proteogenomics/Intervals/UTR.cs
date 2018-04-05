using System.Collections.Generic;

namespace Proteogenomics
{
    public abstract class UTR :
        Interval
    {
        protected UTR(Exon parent, string chromID, string strand, long oneBasedStart, long oneBasedEnd, HashSet<Variant> variants)
            : base(parent, chromID, strand, oneBasedStart, oneBasedEnd, variants)
        {
        }

        protected UTR(UTR utr)
            : base(utr)
        {
        }

        public abstract bool is3Prime();

        public abstract bool is5Prime();

        public abstract override bool CreateVariantEffect(Variant variant, VariantEffects variantEffects);
    }
}