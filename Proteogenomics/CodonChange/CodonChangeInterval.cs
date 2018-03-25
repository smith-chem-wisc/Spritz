namespace Proteogenomics
{
    public class CodonChangeInterval
        : CodonChange
    {
        public CodonChangeInterval(Variant seqChange, Transcript transcript, VariantEffects changeEffects)
            : base(seqChange, transcript, changeEffects)
        {
            ReturnNow = false; // An interval may affect more than one exon
        }

        /// <summary>
        /// Interval is not a variant, nothing to do
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        protected override bool ChangeCodon(Exon exon)
        {
            return false;
        }
    }
}