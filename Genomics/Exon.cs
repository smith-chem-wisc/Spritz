using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class Exon : ChromosomeSegment
    {
       public Exon(string id, Chromosome chrom, string strand, int start, int end, string name, string biotype) 
            : base(id, chrom, strand, start, end, name, biotype)
        { }
    }
}
