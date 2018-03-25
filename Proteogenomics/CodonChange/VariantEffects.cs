using System;
using System.Collections.Generic;
using System.Text;

namespace Proteogenomics
{
    public class VariantEffects
    {
        public List<VariantEffect> Effects { get; set; } = new List<VariantEffect>();

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        /// <param name="effectImpact"></param>
        /// <param name="message"></param>
        public void add(Variant variant, Interval marker, EffectType effectType, EffectImpact effectImpact, String message)
        {
            VariantEffect effNew = new VariantEffect(variant);
            effNew.set(marker, effectType, effectImpact, message);
            add(effNew);
        }

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        /// <param name="message"></param>
        public void add(Variant variant, Interval marker, EffectType effectType, String message)
        {
            add(variant, marker, effectType, VariantEffect.EffectDictionary[effectType], message);
        }

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="variantEffect"></param>
        public void add(VariantEffect variantEffect)
        {
            Effects.Add(variantEffect);
        }

        /// <summary>
        /// Add: If possible, only add an effect type (otherwise add the full effect)
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        public void addEffectType(Variant variant, Interval marker, EffectType effectType)
        {
            if (canAddType(variant, marker))
            {
                get().addEffect(effectType);
            }
            else add(variant, marker, effectType, VariantEffect.EffectDictionary[effectType], "");
        }

        public void addErrorWarning(Variant variant, ErrorWarningType errwarn)
        {
            VariantEffect veff = get();
            if (veff != null) veff.addErrorWarningInfo(errwarn);
            else
            {
                veff = new VariantEffect(variant);
                veff.addErrorWarningInfo(errwarn);
                add(veff);
            }
        }

        /// <summary>
        /// Can we add an effectType to the previous variatnEffect?
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <returns>true if transcript IDs and variant's genotypes match (i.e. we can add effectType)</returns>
        private bool canAddType(Variant variant, Interval marker)
        {
            VariantEffect veff = get();
            if (veff == null || veff.getVariant() == null)
            {
                return false;
            }

            // Do genotypes match?
            var gt = veff.getVariant().Genotype;
            var vgt = variant.Genotype;
            if (((vgt != null) ^ (gt != null)) // One null and one non-null?
                    || ((vgt != null) && (gt != null) && !variant.Genotype.Equals(variant.Genotype)) // Both non-null, but different?
            )
            {
                return false;
            }

            // Do transcripts match?
            Transcript trMarker = (Transcript)marker.findParent(typeof(Transcript));

            Transcript tr = veff.getTranscript();
            if (tr == null || trMarker == null) { return false; }

            return tr.ID == trMarker.ID;
        }

        /// <summary>
        /// Get (or create) the latest ChangeEffect
        /// </summary>
        /// <returns></returns>
        public VariantEffect get()
        {
            if (Effects.Count == 0)
            {
                return null;
            }
            return Effects[Effects.Count - 1];
        }

        public VariantEffect get(int index)
        {
            return Effects[index];
        }

        public bool hasMarker()
        {
            VariantEffect veff = get();
            if (veff == null) return false;
            return veff.getMarker() != null;
        }

        public bool isEmpty()
        {
            return Effects.Count == 0;
        }

        public void setMarker(Interval marker)
        {
            VariantEffect veff = get();
            if (veff != null) { veff.setMarker(marker); }
        }

        public int size()
        {
            return Effects.Count;
        }

        public void sort()
        {
            Effects.Sort();
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append("Effects; " + size() + "\n");

            foreach (VariantEffect eff in Effects)
            {
                sb.Append("\t" + eff.toStr() + "\n");
            }
            return sb.ToString();
        }
    }
}