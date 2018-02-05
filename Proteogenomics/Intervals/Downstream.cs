using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteogenomics
{
    public class Downstream
        : Interval
    {

        #region Constructors

        public Downstream(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {

        }

        public Downstream(Downstream downstream) 
            : base(downstream)
        {

        }

        #endregion Constructors

    }
}
