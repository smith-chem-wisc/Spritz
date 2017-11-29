using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class ExtensionMethods
    {
        public static IEnumerable<IEnumerable<T>> Combinations<T>(this IEnumerable<T> elements, int k)//given an array of elements, it returns all combination sub arrays of length k
        {
            return k == 0 ?
                new[] { new T[0] } :
                elements.SelectMany((e, i) => elements.Skip(i + 1).Combinations(k - 1).Select(c => (new[] { e }).Concat(c)));
        }
    }
}
