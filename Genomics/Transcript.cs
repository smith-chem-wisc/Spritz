using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class Transcript : ChromosomeSegment
    {
        public Gene gene { get; set; }
        public List<SpliceJunction> splices { get; set; } = new List<SpliceJunction>();
        public List<Exon> exons { get; set; } = new List<Exon>();
        public int start_codon_start { get; set; }
        public int stop_codon_start { get; set; }

        public Transcript(string id, Chromosome chrom, string strand, int trans_start, int trans_end, Gene gene, string name, string biotype, int start_codon_start, int stop_codon_start) 
            : base(id, chrom, strand, trans_start, trans_end, name, biotype)
        {
            this.gene = gene;
            this.start_codon_start = start_codon_start;
            this.stop_codon_start = stop_codon_start;
        }

        public Transcript(string id, Chromosome chrom, string strand, int trans_start, int trans_end, Gene gene) 
            : this(id, chrom, strand, trans_start, trans_end, gene, "", "", -1, -1)
        {
            this.gene = gene;
        }

        public void add_exon(Exon exon)
        {
            this.exons.Add(exon);
            if (!(from gene_exon in this.gene.exons select exon.@equals(gene_exon)).Contains(true))
                this.gene.exons.Add(exon);
        }

        public bool valid_exons()
        {
            //List<Exon> exons_by_start = ;
            return false;
        }
    }
}
