namespace Proteogenomics
{
    public class CodonChangeInv
        : CodonChange
    {
        public CodonChangeInv(Variant variant, Transcript transcript, VariantEffects variantEffects)
            : base(variant, transcript, variantEffects)
        {
        }

        public override void ChangeCodon()
        {
            if (Variant.Includes(Transcript))
            {
                // Whole transcript inverted?
                EffectNoCodon(Transcript, EffectType.TRANSCRIPT_INVERSION);
            }
            else
            {
                // Part of the transcript is inverted

                // Does the inversion affect any exon?
                bool intersectsExons = false;
                foreach (Exon ex in Transcript.Exons)
                {
                    if (Variant.Intersects(ex))
                    {
                        intersectsExons = true;
                        break;
                    }
                }

                // Annotate
                if (intersectsExons) { exons(); }
                else { intron(); }
            }
        }

        /// <summary>
        /// One or more exons fully included (no partial overlap)
        /// </summary>
        private void exons()
        {
            Interval cdsMarker = null;
            if (Transcript.IsProteinCoding())
            {
                cdsMarker = Transcript.CdsMarker();
            }

            foreach (Exon ex in Transcript.Exons)
            {
                if (Variant.Intersects(ex))
                {
                    EffectImpact impact = EffectImpact.LOW;

                    // Is the variant affecting a coding part of the exon?
                    // If so, then this is a HIGH impact effect.
                    if (cdsMarker != null && Variant.Intersect(ex).Intersects(cdsMarker))
                    {
                        impact = EffectImpact.HIGH;
                    }

                    // Is the whole exon inverted or just part of it?
                    EffectType effType = Variant.Includes(ex) ? EffectType.EXON_INVERSION : EffectType.EXON_INVERSION_PARTIAL;

                    EffectNoCodon(ex, effType, impact);
                }
            }
        }

        /// <summary>
        /// Inversion does not intersect any exon
        /// </summary>
        private void intron()
        {
            EffectNoCodon(Transcript, EffectType.INTRON);
        }
    }
}