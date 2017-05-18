using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class Transcript : ChromosomeSegment
    {

        #region Public Properties

        public string Name { get; set; }
        public string Biotype { get; set; }
        public Gene Gene { get; set; }
        public List<SpliceJunction> splices { get; set; } = new List<SpliceJunction>();
        public List<Exon> exons { get; set; } = new List<Exon>();
        public int start_codon_start { get; set; } = -1;
        public int stop_codon_start { get; set; } = -1;

        #endregion Public Properties

        #region Public Constructors

        public Transcript(string id, Chromosome chrom, string strand, int trans_start, int trans_end, string biotype, Gene gene)
            : base(id, chrom, strand, trans_start, trans_end)
        {
            this.Biotype = biotype;
            this.Gene = gene;
        }

        #endregion Public Constructors

        #region Public Methods

        public void add_exon(Exon exon)
        {
            exons.Add(exon);
            if (!Gene.exons.Any(x => exon.@equals(x)))
                Gene.exons.Add(exon);
        }

        public bool valid_exons()
        {
            //List<Exon> exons_by_start = ;
            return false;
        }

        #endregion Public Methods

    }
}
