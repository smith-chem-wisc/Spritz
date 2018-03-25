namespace Proteogenomics
{
    public class CodonChangeDup
        : CodonChangeStructural
    {
        public CodonChangeDup(Variant variant, Transcript transcript, VariantEffects variantEffects)
            : base(variant, transcript, variantEffects)
        {
            coding = transcript.isProteinCoding();// || Config.get().isTreatAllAsProteinCoding();
        }

        /**
         * Analyze whether the duplication is past transcript's coding region
         *
         * E.g.:   Transcript  = chr1:100-200
         *         Duplication = chr1:150-999
         *         The duplicated segment starts at base 1000, which is beyond's
         *         transcript's end, so it probably has no effect on the amino
         *         acid sequence
         *
         *  Rationale:
         *  If we have two genes:
         *
         *     | gene_1 |                 |          gene_2            |
         *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
         *        |___________dup____________|
         *
         *  Then this duplication seems to disrupt gene_2:
         *
         *     | gene_1 |                 |          gene_2                                        |
         *  ---<<<<<<<<<<----------------->>>><<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
         *        |___________dup____________||___________dup____________|
         *
         *  Whereas this one does not, because the duplication affects the gene
         *  after the gene's coding region:
         *
         *     | gene_1 |                 |          gene_2            |
         *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>-----
         *                                     |___________dup____________|
         *
         *     | gene_1 |                 |          gene_2            |
         *  ---<<<<<<<<<<----------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>--->>>>>>>>>>>>>>>>>>>>>>>>>-----
         *                                     |___________dup____________||___________dup____________|
         *
         * @return true if the duplication is beyond transcript's end
         */

        private bool BeyondTranscript()
        {
            if (coding)
            {
                if (Transcript.Strand == "+") return Variant.OneBasedEnd > Transcript.cdsEnd;
                return Variant.OneBasedEnd > Transcript.cdsStart;
            }

            return Variant.OneBasedEnd > Transcript.OneBasedEnd;
        }

        protected override void EffectTranscript()
        {
            EffectNoCodon(Transcript, EffectType.TRANSCRIPT_DUPLICATION);
        }

        protected override void Exons()
        {
            if (BeyondTranscript())
            {
                // Is the effect of a duplication beyond transcript's end?
                // Then it probably does not have much impact
                EffectImpact impact = coding ? EffectImpact.LOW : EffectImpact.MODIFIER;
                if (exonFull > 0) EffectNoCodon(Transcript, EffectType.EXON_DUPLICATION, impact);
                if (exonPartial > 0) EffectNoCodon(Transcript, EffectType.EXON_DUPLICATION_PARTIAL, impact);
                return;
            }

            if (coding) ExonsCoding();
            else ExonsNoncoding();
        }

        /// <summary>
        /// One or more exons fully included (no partial overlap)
        /// </summary>
        protected override void ExonsCoding()
        {
            CodonsRefAlt();

            if (exonFull > 0) Effect(Transcript, EffectType.EXON_DUPLICATION, false);
            if (exonPartial > 0) Effect(Transcript, EffectType.EXON_DUPLICATION_PARTIAL, false);

            // Is this duplication creating a frame-shift?
            int lenDiff = cdsAlt.Length - cdsRef.Length;
            if (lenDiff % 3 != 0) Effect(Transcript, EffectType.FRAME_SHIFT, false);
        }

        /// <summary>
        /// Effects for non-coding transcripts
        /// </summary>
        protected override void ExonsNoncoding()
        {
            if (exonFull > 0) EffectNoCodon(Transcript, EffectType.EXON_DUPLICATION, EffectImpact.MODIFIER);
            if (exonPartial > 0) EffectNoCodon(Transcript, EffectType.EXON_DUPLICATION_PARTIAL, EffectImpact.MODIFIER);
        }

        /// <summary>
        /// Inversion does not intersect any exon
        /// </summary>
        protected override void Intron()
        {
            EffectNoCodon(Transcript, EffectType.INTRON);
        }
    }
}