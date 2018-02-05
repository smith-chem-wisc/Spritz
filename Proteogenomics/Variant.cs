using Bio.VCF;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Variant
        : VariantContext
    {

        #region Public Properties

        /// <summary>
        /// The (start) position of this variation
        /// </summary>
        public int OneBasedStart { get; set; }

        /// <summary>
        /// The reference allele
        /// </summary>
        public Allele ReferenceAllele { get; set; }

        /// <summary>
        /// The alternate allele
        /// </summary>
        public Allele AlternateAllele { get; set; }

        /// <summary>
        /// The reference allele as a string
        /// </summary>
        public string ReferenceAlleleString { get; set; }

        /// <summary>
        /// The alternate allele as a string
        /// </summary>
        public string AlternateAlleleString { get; set; }

        /// <summary>
        /// The genotype of this variation
        /// </summary>
        public Genotype Genotype { get; set; }

        /// <summary>
        /// The type of variation: heterozygous/homozygous
        /// </summary>
        public GenotypeType GenotypeType { get; set; }

        /// <summary>
        /// Depth across this position
        /// </summary>
        public int ReadDepth { get; set; }

        /// <summary>
        /// Depth of reads across this position containing the reference allele
        /// </summary>
        public int ReferenceAlleleDepth { get; set; }

        /// <summary>
        /// Depth of reads across this position containing the alternate allele
        /// </summary>
        public int AlternateAlleleDepth { get; set; }

        /// <summary>
        /// Text annotation for protein entries
        /// </summary>
        public string Annotation { get; set; }

        /// <summary>
        /// SnpEff annotaitons for this location.
        /// </summary>
        public List<SnpEffAnnotation> SnpEffAnnotations { get; } = new List<SnpEffAnnotation>();

        #endregion Public Properties

        #region Public Constructor

        public Variant(VariantContext variant, int zeroBasedAlleleIndex)
            : base(variant)
        {
            OneBasedStart = variant.Start;
            ReferenceAllele = variant.Reference;
            ReferenceAlleleString = variant.Reference.BaseString;
            Genotype = variant.Genotypes.First(); // assumes just one sample
            GenotypeType = Genotype.Type;
            ReadDepth = Genotype.DP;
            ReferenceAlleleDepth = Genotype.AD[0];
            AlternateAlleleDepth = Genotype.AD[1];
            AlternateAllele = Genotype.GetAllele(1);
            AlternateAlleleString = variant.AlternateAlleles[zeroBasedAlleleIndex].BaseString;

            if (this.Attributes.TryGetValue("ANN", out object snpEffAnnotationObj))
            {
                if (snpEffAnnotationObj as string[] != null) SnpEffAnnotations = (snpEffAnnotationObj as string[]).Select(x => new SnpEffAnnotation(this, x)).ToList();
                if (snpEffAnnotationObj as string != null) SnpEffAnnotations = new List<SnpEffAnnotation> { new SnpEffAnnotation(this, snpEffAnnotationObj as string) };
            }
        }

        #endregion Public Constructor

    }
}
