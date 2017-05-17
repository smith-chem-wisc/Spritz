using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Proteomics;

namespace Genomics
{
    public class Chromosome
    {
        public string name { get; set; }
        public string sequence { get; set; }
        public List<Gene> genes { get; set; }
        public List<AminoAcidPolymer> amino_acid_sequences { get; set; }
        public int length { get { return sequence.Length; } }

        #region Public Constructor

        public Chromosome(string name)
        {
            this.name = name;
        } 

        #endregion

        #region Public Methods

        public bool contains(string gene_name)
        {
            return genes.Select(g => g.name).Contains(gene_name);
        }

        public void sort()
        {
            //TODO
        }

        public Gene get_gene_by_position(int position)
        {
            return new Gene("", new Chromosome(""), "", 0, 0);
        }

        public override string ToString()
        {
            return this.sequence;
        }

        #endregion Public Methods
    }
}
