using Genomics;
using NUnit.Framework;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class IUPACTesting
    {
        [Test]
        public void ListsContainOnlyCanonicalMonomers()
        {
            char[] aa_3to1_values = new char[NucleotideSequence.amino_acids_3to1.Values.Count];
            NucleotideSequence.amino_acids_3to1.Values.CopyTo(aa_3to1_values, 0);
            foreach (char c in NucleotideSequence.amino_acids)
            {
                //Assert.IsTrue(Individual.amino_acids_1to3.TryGetValue(c, null));
                //Assert.IsTrue(aa_3to1_values.)
            }
        }
    }
}
