namespace Proteogenomics
{
    public class UTR3Prime
        : UTR
    {
        #region Constructors

        public UTR3Prime(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public UTR3Prime(UTR3Prime utr)
            : base(utr)
        {
        }

        #endregion Constructors

        #region Public Methods

        public override bool is3Prime()
        {
            return true;
        }

        public override bool is5Prime()
        {
            return false;
        }

        #endregion Public Methods
    }
}