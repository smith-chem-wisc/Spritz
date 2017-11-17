using Bio;
using Bio.IO.FastA;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Genome
    {

        #region Public Properties

        public List<ISequence> Chromosomes { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public Genome(string genomeFastaLocation)
        {
            Chromosomes = new FastAParser().Parse(genomeFastaLocation).ToList();
        }

        #endregion Public Constructor

    }
}
