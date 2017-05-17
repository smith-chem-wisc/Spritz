using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class Insertion : Indel
    {

        #region Public Constructor

        public Insertion(Chromosome chrom, int position, string id, string reference, string alternate, double qual, string filter, Dictionary<string, string> info)
            : base(chrom, position, id, reference, alternate, qual, filter, info)
        { }

        #endregion Public Constructor

    }
}
