using Bio.VCF;
using System.Collections.Generic;

namespace Genomics
{
    public class VariantLite
    {
        public string ReferenceAllele { get; set; }
        public string AlternateAllele { get; set; }
        public double AlleleFrequency { get; set; }
        public int OneBasedStart { get; set; }

        public VariantLite(VariantContext variant, int zeroBasedAlleleIndex)
        {
            OneBasedStart = variant.Start;
            ReferenceAllele = variant.Reference.BaseString;
            AlternateAllele = variant.AlternateAlleles[zeroBasedAlleleIndex].BaseString;
            if (double.TryParse(variant.GetAttributeAsString("AF", "").Split(',')[zeroBasedAlleleIndex], out double af)) AlleleFrequency = af;
            else AlleleFrequency = 1;
        }

        public static IEnumerable<VariantLite> ParseVariantContext(VariantContext variant)
        {
            for (int i = 0; i < variant.AlternateAlleles.Count; i++)
            {
                yield return new VariantLite(variant, i);
            }
        }
    }
}
