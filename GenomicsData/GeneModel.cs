using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Genomics;

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
