using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class Chromosome
    {

        #region Public Constructors

        public string Name { get; set; }
        public string Sequence { get; set; }
        public int Length { get { return Sequence.Length; } }

        #endregion Public Constructors

        #region Public Constructor

        public Chromosome(string name, string sequence)
        {
            this.Name = name;
            this.Sequence = sequence;
        } 

        #endregion

    }
}
