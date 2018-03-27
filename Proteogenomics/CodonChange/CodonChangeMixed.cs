//using System;
//using System.Collections.Generic;
//using System.Linq;

//namespace Proteogenomics
//{
//    public class CodonChangeMixed
//        : CodonChangeMnv
//    {
//        public static bool debug = false;

//        private int oldCodonCdsStart = -1;
//        private int oldCodonCdsEnd = -1;

//        private Variant mnp;
//        private Variant indel;
//        private CodonChange codonChangeMnp;
//        private CodonChange codonChangeIndel;
//        private VariantEffects variantEffectsOri;

//        public CodonChangeMixed(Variant variant, Transcript transcript, VariantEffects variantEffects)
//            : base(variant, transcript, variantEffects)
//        {
//            ReturnNow = false;
//            RequireNetCdsChange = false;

//            // Decompose mixed variant into 'SNP/MNP' + 'INS/DEL'
//            int minLen = Math.Min(variant.ReferenceAlleleString.Length, variant.AlternateAlleleString.Length);

//            string reff = variant.ReferenceAlleleString;
//            string refMnp = reff.Substring(0, minLen);
//            string refIndel = reff.Substring(minLen);

//            string alt = variant.AlternateAlleleString;
//            string altMnp = alt.Substring(0, minLen);
//            string altIndel = alt.Substring(minLen);

//            mnp = new Variant(variant.VariantContext, variant.Chromosome, variant.VariantContext.AlternateAlleles.IndexOf(variant.VariantContext.AlternateAlleles.FirstOrDefault(a => a.BaseString == variant.AlternateAlleleString)));
//            indel = new Variant(variant.Chromosome, variant.OneBasedStart + minLen, refIndel, altIndel, variant.getId());

//            // Create codon changes
//            variantEffectsOri = variantEffects;
//            VariantEffects = new VariantEffects();
//            codonChangeMnp = Factory(mnp, transcript, this.VariantEffects);
//            codonChangeIndel = Factory(indel, transcript, this.VariantEffects);
//        }

//        public override void ChangeCodon()
//        {
//            codonOldNew();

//            codonChangeMnp.ChangeCodon();
//            codonChangeIndel.ChangeCodon();

//            // Set highest impact variant effect
//            if (VariantEffects.Effects.Count == 0 || VariantEffects.Effects.Count == 0) return; // Nothing to do

//            VariantEffects.Effects.Sort();
//            VariantEffect varEff = VariantEffects.Effects[0];

//            // Add main effect
//            varEff = Effect(varEff.getMarker(), varEff.getEffectType(), false);

//            // Add 'additional' effects
//            for (int i = 0; i < VariantEffects.Effects.Count; i++)
//            {
//                List<EffectType> effTypes = VariantEffects.Effects[i].effectTypes;
//                for (int j = 0; j < effTypes.Count; j++)
//                {
//                    EffectType effType = effTypes[j];
//                    if (!varEff.hasEffectType(effType)) varEff.addEffect(effType);
//                }
//            }

//            variantEffectsOri.add(varEff);
//        }

//        private void codonNum()
//        {
//            if (Transcript.IsStrandPlus())
//            {
//                CodonStartNumber = codonChangeMnp.CodonStartNumber;
//                CodonStartIndex = codonChangeMnp.CodonStartIndex;
//            }
//            else
//            {
//                CodonStartNumber = codonChangeIndel.CodonStartNumber;
//                CodonStartIndex = codonChangeIndel.CodonStartIndex;
//            }
//        }
//    }
//}