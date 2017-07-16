using System.Collections.Generic;

namespace Genomics
{
    public class Sample
    {

        #region Public Properties

        public string name { get; set; }

        public Dictionary<string, Chromosome> chroms { get; set; }

        public List<SequenceVariant> sequence_variants { get; set; }

        public List<LocalHaplotype> local_haplotypes { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public Sample(string name)
        {
            this.name = name;
            chroms = new Dictionary<string, Chromosome>();
            sequence_variants = new List<SequenceVariant>();
            local_haplotypes = new List<LocalHaplotype>();
        }

        #endregion Public Constructor

    }
}
