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

        public override bool CreateVariantEffect(Variant variant, VariantEffects variantEffects)
        {
            if (!Intersects(variant)) return false;

            //for (SpliceSite ss : spliceSites)
            //    if (ss.intersects(variant)) ss.variantEffect(variant, variantEffects);

            // Add intron part
            variantEffects.AddEffectType(variant, this, EffectType.INTRON);

            return true;
        }
    }
}