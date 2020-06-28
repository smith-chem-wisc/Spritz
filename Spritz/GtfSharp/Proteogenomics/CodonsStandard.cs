using System.Collections.Generic;

namespace Proteogenomics
{
    public static class CodonsStandard
    {
        /// <summary>
        /// All start codons are translated as "M", even if alternative start codons, since an different tRNA is used. https://en.wikipedia.org/wiki/Start_codon
        /// </summary>
        public static readonly string DEFAULT_START_AA = "M";

        /// <summary>
        /// Start codons for the standard genetic code.
        /// </summary>
        public static readonly HashSet<string> START_CODONS = new HashSet<string> { "ATG", "TTG", "CTG" };
    }
}