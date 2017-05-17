using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class Gene : ChromosomeSegment
    {
        public List<Transcript> transcripts { get; set; } = new List<Transcript>();
        public List<SpliceJunction> splices { get; set; } = new List<SpliceJunction>();
        public HashSet<Exon> exons { get; set; } = new HashSet<Exon>();
        public Gene(string id, Chromosome chrom, string strand, int gene_start, int gene_end, string name, string biotype) 
            : base(id, chrom, strand, gene_start, gene_end, name, biotype)
        { }
        public Gene(string id, Chromosome chrom, string strand, int gene_start, int gene_end) 
            : this(id, chrom, strand, gene_start, gene_end, "", "")
        { }

        public Transcript get_fully_spliced_transcript()
        {
            Transcript fully_spliced = new Transcript(this.name, this.chrom, this.strand, this.start, this.end, this);
            fully_spliced.exons = exons.ToList();
            fully_spliced.valid_exons();
            return fully_spliced;
        }
    }
}
