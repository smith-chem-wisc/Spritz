using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
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
        /// <param name="message"></param>
        public void AddEffect(Variant variant, Interval marker, EffectType effectType, String message)
        {
            AddEffect(variant, marker, effectType, EffectTypeMethods.EffectDictionary[effectType], message);
        }

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="variantEffect"></param>
        public void AddEffect(VariantEffect variantEffect)
        {
            Effects.Add(variantEffect);
        }

        /// <summary>
        /// Add: If possible, only add an effect type (otherwise add the full effect)
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        public void AddEffectType(Variant variant, Interval marker, EffectType effectType)
        {
            if (CanAddType(variant, marker))
            {
                Get().AddEffect(effectType);
            }
            else
            {
                AddEffect(variant, marker, effectType, EffectTypeMethods.EffectDictionary[effectType], "");
            }
        }

        /// <summary>
        /// Add an error or warning
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="errwarn"></param>
        public void AddErrorWarning(Variant variant, ErrorWarningType errwarn)
        {
            VariantEffect veff = Get();
            if (veff != null)
            {
                veff.AddErrorWarningInfo(errwarn);
            }
            else
            {
                veff = new VariantEffect(variant);
                veff.AddErrorWarningInfo(errwarn);
                AddEffect(veff);
            }
        }

        /// <summary>
        /// Get (or create) the latest ChangeEffect
        /// </summary>
        /// <returns></returns>
        public VariantEffect Get()
        {
            return Effects.Count == 0 ? null : Effects[Effects.Count - 1];
        }

        public bool HasMarker()
        {
            VariantEffect veff = Get();
            if (veff == null) return false;
            return veff.Marker != null;
        }

        public void SetMarker(Interval marker)
        {
            VariantEffect veff = Get();
            if (veff != null) veff.SetMarker(marker);
        }

        /// <summary>
        /// Get a string representing the annotation for this variant effect on a transript
        /// </summary>
        /// <returns></returns>
        public string TranscriptAnnotation()
        {
            StringBuilder sb = new StringBuilder();
            Variant theVariant = Effects[0].Variant;
            sb.Append("variant:" + theVariant.ToString() + " ");
            foreach (VariantEffect eff in Effects)
            {
                sb.Append(eff.TranscriptAnnotation());
            }
            return sb.ToString();
        }

        /// <summary>
        /// Gets a list of the protein sequence variations noted as effects (has length of one in reference to a single transcript)
        /// </summary>
        /// <returns></returns>
        public List<SequenceVariation> ProteinSequenceVariation()
        {
            return Effects.Where(eff => eff.IsNonsynonymous()) // has coding effect
                .Select(eff => new SequenceVariation(eff.CodonNum + 1, eff.CodonNum + eff.AlternateAA.Length, eff.ReferenceAA, eff.AlternateAA, eff.TranscriptAnnotation())).ToList();
        }

        /// <summary>
        /// Get a string representing this object
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append("Effects; " + Effects.Count.ToString() + "\n");

            foreach (VariantEffect eff in Effects)
            {
                sb.Append("\t" + eff.ToStr() + "\n");
            }
            return sb.ToString();
        }

        /// <summary>
        /// Add an effect
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <param name="effectType"></param>
        /// <param name="effectImpact"></param>
        /// <param name="message"></param>
        private void AddEffect(Variant variant, Interval marker, EffectType effectType, EffectImpact effectImpact, String message)
        {
            VariantEffect effNew = new VariantEffect(variant);
            effNew.Set(marker, effectType, effectImpact, message);
            AddEffect(effNew);
        }

        /// <summary>
        /// Can we add an effectType to the previous variatnEffect?
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="marker"></param>
        /// <returns>true if transcript IDs and variant's genotypes match (i.e. we can add effectType)</returns>
        private bool CanAddType(Variant variant, Interval marker)
        {
            VariantEffect veff = Get();
            if (veff == null || veff.Variant == null)
            {
                return false;
            }

            // Do genotypes match?
            var gt = veff.Variant.Genotype;
            var vgt = variant.Genotype;
            if (((vgt != null) ^ (gt != null)) // One null and one non-null?
                    || ((vgt != null) && (gt != null) && !variant.Genotype.Equals(variant.Genotype)) // Both non-null, but different?
            )
            {
                return false;
            }

            // Do transcripts match?
            Transcript trMarker = (Transcript)marker.FindParent(typeof(Transcript));

            Transcript tr = veff.GetTranscript();
            if (tr == null || trMarker == null) { return false; }

            return tr.ID == trMarker.ID;
        }
    }
}