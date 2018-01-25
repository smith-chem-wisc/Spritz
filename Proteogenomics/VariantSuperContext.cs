using Bio.VCF;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class VariantSuperContext 
        : VariantContext
    {

        #region Public Properties

        /// <summary>
        /// Custom variant object, which has one allele each
        /// </summary>
        public List<Variant> Variants { get; } = new List<Variant>();

        /// <summary>
        /// SnpEff annotaitons for this location.
        /// </summary>
        public List<SnpEffAnnotation> SnpEffAnnotations { get; } = new List<SnpEffAnnotation>();

        #endregion Public Properties

        #region Public Constructor

        public VariantSuperContext(VariantContext variantContext) 
            : base(variantContext)
        {
            Variants = Proteogenomics.Variant.ParseVariantContext(variantContext).ToList();

            if (this.Attributes.TryGetValue("ANN", out object snpEffAnnotationObj))
            {
                if (snpEffAnnotationObj as string[] != null) SnpEffAnnotations = (snpEffAnnotationObj as string[]).Select(x => new SnpEffAnnotation(this, x)).ToList();
                if (snpEffAnnotationObj as string != null) SnpEffAnnotations = new List<SnpEffAnnotation> { new SnpEffAnnotation(this, snpEffAnnotationObj as string) };
            }
        }

        #endregion Public Constructor

    }
}
