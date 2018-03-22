using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class CodonChangeMixed
        : CodonChangeMnv
    {
        public static bool debug = false;

        private int oldCodonCdsStart = -1;
        private int oldCodonCdsEnd = -1;

        private Variant mnp;
        private Variant indel;
        private CodonChange codonChangeMnp;
        private CodonChange codonChangeIndel;
        private List<VariantEffect> variantEffectsOri;

        public CodonChangeMixed(Variant variant, Transcript transcript, List<VariantEffect> variantEffects)
            : base(variant, transcript, variantEffects)
        {
            ReturnNow = false;
            RequireNetCdsChange = false;

            // Decompose mixed variant into 'SNP/MNP' + 'INS/DEL'
            int minLen = Math.Min(variant.getReference().length(), variant.getAlt().length());

            string reff = variant.getReference();
            string refMnp = reff.Substring(0, minLen);
            string refIndel = reff.Substring(minLen);

            string alt = variant.getAlt();
            string altMnp = alt.Substring(0, minLen);
            string altIndel = alt.Substring(minLen);

            mnp = new Variant(variant.getChromosome(), variant.getStart(), refMnp, altMnp, variant.getId());
            indel = new Variant(variant.getChromosome(), variant.getStart() + minLen, refIndel, altIndel, variant.getId());

            // Create codon changes
            variantEffectsOri = variantEffects;
            this.VariantEffects = new List<VariantEffect>();
            codonChangeMnp = Factory(mnp, transcript, this.VariantEffects);
            codonChangeIndel = Factory(indel, transcript, this.VariantEffects);
        }

        public override void ChangeCodon()
        {
            codonOldNew();

            codonChangeMnp.ChangeCodon();
            codonChangeIndel.ChangeCodon();

            // Set highest impact variant effect
            if (VariantEffects.Count == 0 || VariantEffects.Count == 0) return; // Nothing to do

            VariantEffects.Sort();
            VariantEffect varEff = VariantEffects[0];

            // Add main effect
            varEff = Effect(varEff.getMarker(), varEff.getEffectType(), false);

            // Add 'additional' effects
            for (int i = 0; i < VariantEffects.Count; i++)
            {
                List<EffectType> effTypes = VariantEffects[i].getEffectTypes();
                for (int j = 0; j < effTypes.Count; j++)
                {
                    EffectType effType = effTypes[j];
                    if (!varEff.hasEffectType(effType)) varEff.addEffect(effType);
                }
            }

            variantEffectsOri.Add(varEff);
        }

        private void codonNum()
        {
            if (Transcript.Strand == "+")
            {
                CodonStartNumber = codonChangeMnp.CodonStartNumber;
                CodonStartIndex = codonChangeMnp.CodonStartIndex;
            }
            else
            {
                CodonStartNumber = codonChangeIndel.CodonStartNumber;
                CodonStartIndex = codonChangeIndel.CodonStartIndex;
            }
        }
    }
}