using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace Proteogenomics
{
    public abstract class UTR :
        Interval
    {

        #region Constructors

        public UTR(string chromID, string strand, long oneBasedStart, long oneBasedEnd) 
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {

        }

        public UTR(UTR utr) 
            : base(utr)
        {

        }

        #endregion Constructors

        #region Public Methods

        public abstract bool is3Prime();

        public abstract bool is5Prime();

        #endregion Public Methods

    }
}
