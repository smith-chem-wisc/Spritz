namespace Proteogenomics
{
    public class Intergenic :
        Interval
    {
        #region Public Properties

        public Gene LeftGene { get; set; }

        public Gene RightGene { get; set; }

        #endregion Public Properties

        #region Constructors

        public Intergenic(string chromID, string strand, long oneBasedStart, long oneBasedEnd) :
            base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public Intergenic(Intergenic intergenic) : base(intergenic)
        {
        }

        #endregion Constructors
    }
}