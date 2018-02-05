using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteogenomics
{
    public class Intron
        : Interval
    {

        #region Constructors

        public Intron(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {

        }

        public Intron(Intron intron)
            : base(intron)
        {

        }

        #endregion Constructors

    }
}
