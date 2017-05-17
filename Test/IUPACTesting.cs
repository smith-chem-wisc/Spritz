using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Genomics;

namespace Test
{
    [TestFixture]
    public class IUPACTesting
    {
        [Test]
        public void ListsContainOnlyCanonicalMonomers()
        {
            char[] aa_3to1_values = new char[Sample.amino_acids_3to1.Values.Count];
            Sample.amino_acids_3to1.Values.CopyTo(aa_3to1_values, 0);
            foreach (char c in Sample.amino_acids)
            {
                //Assert.IsTrue(Individual.amino_acids_1to3.TryGetValue(c, null));
                //Assert.IsTrue(aa_3to1_values.)
            }
        }
    }
}
