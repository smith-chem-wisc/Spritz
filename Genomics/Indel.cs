using System.Collections.Generic;

namespace Genomics
{
    public class Indel : SequenceVariant
    {

        #region Public Constructor

        public Indel(Chromosome chrom, int position, string id, string reference, string alternate, double qual, string filter, Dictionary<string, string> info)
            : base(chrom, position, id, reference, alternate, qual, filter, info)
        { }

        #endregion Public Constructor

    }
}
