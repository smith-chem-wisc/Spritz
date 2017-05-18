namespace Genomics
{
    public class Exon : ChromosomeSegment
    {

        #region Public Properties

        Transcript transcript { get; set; }
        Gene gene { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public Exon(string id, Chromosome chrom, string strand, int start, int end, Transcript transcript, Gene gene)
             : base(id, chrom, strand, start, end)
        {
            this.transcript = transcript;
            this.gene = gene;
        }

        #endregion Public Constructor

    }
}
