using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class Gene : ChromosomeSegment
    {

        #region Public Properties

        public string Name { get; set; }
        public string Biotype { get; set; }
        public List<Transcript> transcripts { get; set; } = new List<Transcript>();
        public List<SpliceJunction> splices { get; set; } = new List<SpliceJunction>();
        public HashSet<Exon> exons { get; set; } = new HashSet<Exon>();

        #endregion

        #region Public Constructors

        public Gene(string id, Chromosome chrom, string strand, int gene_start, int gene_end, string name, string biotype)
            : base(id, chrom, strand, gene_start, gene_end)
        { }

        #endregion Public Constructors

        #region Public Methods

        public Transcript get_fully_spliced_transcript()
        {
            Transcript fully_spliced = new Transcript(Name + "_full", Chrom, Strand, OneBasedStart, OneBasedEnd, Biotype, this);
            fully_spliced.exons = exons.ToList();
            fully_spliced.valid_exons();
            return fully_spliced;
        }

        #endregion Public Methods

    }
}
