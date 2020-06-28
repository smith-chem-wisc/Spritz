using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class ProteogenomicsUtility
    {
        public static readonly string EnsemblFastaHeaderDelimeter = " ";

        public static char[] amino_acids = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();

        /// <summary>
        /// IUPAC Codes Not Contained in DotNetBio
        /// </summary>
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
            { 'Y', "Tyr" },

            { '*', "Ter" },
            { '-', "Gap" },
            { 'U', "Sec" },
            { 'O', "Pyl" },

            { 'X', "Xaa" },
            { 'Z', "Glx" },
            { 'B', "Asx" },
            { 'J', "Xle" },
        };

        public static Dictionary<string, char> amino_acids_3to1 =
            amino_acids_1to3.ToList().ToDictionary(
                one_to_three => one_to_three.Value,
                one_to_three => one_to_three.Key);

        public static IEnumerable<IEnumerable<T>> Combinations<T>(this IEnumerable<T> elements, int k)//given an array of elements, it returns all combination sub arrays of length k
        {
            return k == 0 ?
                new[] { new T[0] } :
                elements.SelectMany((e, i) => elements.Skip(i + 1).Combinations(k - 1).Select(c => (new[] { e }).Concat(c)));
        }
    }
}