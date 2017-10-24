using Bio;
using Bio.IO.FastA;
using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class Genome
    {

        #region Public Properties

        public List<ISequence> chroms { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public Genome(string genome_fasta_location)
        {
            chroms = new FastAParser().Parse(genome_fasta_location).ToList();
        }

        #endregion Public Constructor

    }
}
