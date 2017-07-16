using Genomics;
using System.Collections.Generic;

namespace GenomicsData
{
    public class GeneModel
    {

        #region Public Properties

        public Dictionary<string, Chromosome> chromosomes { get; set; }
        public List<Gene> genes { get; set; }

        #endregion Public Properties

        #region Public Constructors

        public GeneModel(Dictionary<string, Chromosome> chromosomes, List<Gene> genes)
        {
            this.chromosomes = chromosomes;
            this.genes = genes;
        }

        #endregion Public Constructors

    }
}
