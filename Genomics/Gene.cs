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
        public HashSet<Exon> exons { get; set; } = new HashSet<Exon>();

        #endregion

        #region Public Constructors

        public Gene(string id, Chromosome chrom, string strand, int gene_start, int gene_end, string name, string biotype)
            : base(id, chrom, strand, gene_start, gene_end)
        { }

        #endregion Public Constructors

        #region Public Methods

        // Kind of stupid because the fully-spliced transcript should be annotated in the GTF file
        //public Transcript get_fully_spliced_transcript()
        //{
        //    Transcript fully_spliced = new Transcript(Name + "_full", Chrom, Strand, OneBasedStart, OneBasedEnd, Biotype, this);
        //    fully_spliced.exons = exons.ToList();
        //    List<Exon> combined = new List<Exon>();
        //    List<Exon> remaining = new List<Exon>();
        //    while (remaining.Count > 0)
        //    {
        //        combined.Add(remaining[0]);
        //        fully_spliced.exons.Add(remaining[0].combine_overlapping_sequences(remaining, out IEnumerable<ChromosomeSegment> ));

        //    }
        //    foreach (Exon x in fully_spliced.exons)
        //    {
        //        x.combine_overlapping_sequences(fully_spliced.exons.Where(xy => x.overlaps(xy)));
        //    }
        //    return fully_spliced;
        //}

        #endregion Public Methods

    }
}
