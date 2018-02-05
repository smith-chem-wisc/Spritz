using Bio.VCF;
using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class VariantOld
        : VariantContext
    {

        #region Public Properties

        public VariantContext OriginalContext { get; set; }

        public string ReferenceAllele { get; set; }

        public string AlternateAllele { get; set; }

        public double AlleleFrequency { get; set; }

        //public int AlleleReadDepth { get; set; }

        public int OneBasedStart { get; set; }

        public string Annotation { get; set; }

        public bool Synonymous { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public VariantOld(VariantContext variant, int zeroBasedAlleleIndex)
            : base(variant)
        {
            OriginalContext = variant;
            OneBasedStart = variant.Start;
            ReferenceAllele = variant.Reference.BaseString;
            AlternateAllele = variant.AlternateAlleles[zeroBasedAlleleIndex].BaseString;
            AlleleFrequency = 1;

            if (variant.Attributes.TryGetValue("AF", out object alleleFrequencyObj))
            {
                if (alleleFrequencyObj as string[] != null && double.TryParse((alleleFrequencyObj as string[])[zeroBasedAlleleIndex], out double af)) AlleleFrequency = af;
                else if (double.TryParse(alleleFrequencyObj as string, out double af2)) AlleleFrequency = af2;
                else AlleleFrequency = 1;
            }

            // untested
            //if (variant.Attributes.TryGetValue("DP", out object readDepthObj))
            //{
            //    if (readDepthObj as string[] != null && double.TryParse((readDepthObj as string[])[zeroBasedAlleleIndex], out double dp)) AlleleReadDepth = dp;
            //    else if (double.TryParse(readDepthObj as string, out double dp2)) AlleleReadDepth = dp2;
            //    else AlleleReadDepth = 1;
            //}
        }

        #endregion Public Constructor

        #region Public Method

        public static IEnumerable<VariantOld> ParseVariantContext(VariantContext variant)
        {
            for (int i = 0; i < variant.AlternateAlleles.Count; i++)
            {
                yield return new VariantOld(variant, i);
            }
        }

        #endregion Public Method

    }
}
