using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
//# A local haplotype is a set of alleles known or predicted to exist on the same chromosome. Since there are two of each
//# chromosome in a haploid genome, certain alleles exist 'adjacent' to one another in a genotype.
//# The term for assigning a local haplotype in variant calling workflows is 'phasing.' The result of phasing is a list
//# of alleles that are in phase or out of phase with one another. When in phase, alleles are expected to exist on the
//# same chromosome. When out of phase, there is no association between the alleles.
//#
//# Accordingly, each LocalHaplotype object contains two lists of SequenceVariants, corresponding to two of the same
//# chromosome in a diploid genome. Since homozygous alleles are always annotated as in phase with the previous
//# heterozygous allele in phasing tools, there should be no issue with redundnancy of haplotypes in a gene model.
//#
//# Multiple sequence variants can exist at a locus, i.e. specified in the same line of a variant call file,
//# so be sure to check for that later.
    public class LocalHaplotype : ChromosomeSegment
    {
        public List<SequenceVariant> ploid1 { get; set; }
        public List<SequenceVariant> ploid2 { get; set; }
        public LocalHaplotype(Chromosome chrom, int start) : base(null, chrom, "+", start, start, null, null)
        { }

        public void add(SequenceVariant ploid1_seqvar, SequenceVariant ploid2_seqvar)
        {
            if (ploid1_seqvar != null) this.ploid1.Add(ploid1_seqvar);
            if (ploid2_seqvar != null) this.ploid2.Add(ploid2_seqvar);
        }

        public void update_end(int end)
        {
            this.end = end;
        }
    }
}
