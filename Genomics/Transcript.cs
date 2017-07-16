using Proteomics;
using System.Collections.Generic;

namespace Genomics
{
    public class Transcript : ChromosomeSegment
    {

        #region Public Properties

        public string Biotype { get; set; }
        public Gene Gene { get; set; }
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

        public Protein translate()
        {
            string seq = "";
            foreach (Exon x in exons)
            {
                int start1base = x.OneBasedStart > start_codon_start ? x.OneBasedStart : start_codon_start;
                int end1base = x.OneBasedEnd < stop_codon_start ? x.OneBasedEnd : stop_codon_start;
                if (start1base > end1base)
                {
                    continue;
                }
                seq += x.Chrom.Sequence.get_range(start1base, end1base);
            }
            NucleotideSequence nt_seq = new NucleotideSequence(seq);
            return nt_seq.get_peptide(Chrom.Name.Contains("M") ? NucleotideSequence.mitochon_code : NucleotideSequence.standard_code, ID);
        }

        #endregion Public Methods

    }
}
