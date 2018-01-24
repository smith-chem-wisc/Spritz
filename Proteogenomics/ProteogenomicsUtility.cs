using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class ProteogenomicsUtility
    {

        #region IUPAC Codes Not Contained in DotNetBio

        public static char[] amino_acids = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();

        public static Dictionary<char, string> amino_acids_1to3 = new Dictionary<char, string>()
        {
            { 'A', "Ala" },
            { 'C', "Cys" },
            { 'D', "Asp" },
            { 'E', "Glu" },
            { 'F', "Phe" },
            { 'G', "Gly" },
            { 'H', "His" },
            { 'I', "Ile" },
            { 'K', "Lys" },
            { 'L', "Leu" },
            { 'M', "Met" },
            { 'N', "Asn" },
            { 'P', "Pro" },
            { 'Q', "Gln" },
            { 'R', "Arg" },
            { 'S', "Ser" },
            { 'T', "Thr" },
            { 'V', "Val" },
            { 'W', "Trp" },
            { 'Y', "Tyr" }
        };

        public static Dictionary<string, char> amino_acids_3to1 =
            amino_acids_1to3.ToList().ToDictionary(
                one_to_three => one_to_three.Value,
                one_to_three => one_to_three.Key);

        #endregion IUPAC Codes Not Contained in DotNetBio

    }
}
