namespace Proteogenomics
{
    public class UTR5Prime
        : UTR
    {
        #region Constructors

        public UTR5Prime(string chromID, string strand, long oneBasedStart, long oneBasedEnd) :
            base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public UTR5Prime(UTR5Prime utr) : base(utr)
        {
        }

        #endregion Constructors

        #region Public Methods

        public override bool is3Prime()
        {
            return false;
        }

        public override bool is5Prime()
        {
            return true;
        }

        #endregion Public Methods
    }
}