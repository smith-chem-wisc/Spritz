using NUnit.Framework;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace Test
{
    [TestFixture]
    public class IndependentWrapperTests
    {
        [Test]
        public void StrandAmbiguityFiltering()
        {
            var stringtie = new StringtieWrapper();
            string f = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf_strandProblem.gtf");
            string filtered = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf_strandProblem.filtered.gtf");
            stringtie.FilterGtfEntriesWithoutStrand(f, filtered, false);

            // assert each one has a strand
            var lines = File.ReadAllLines(filtered);
            Assert.IsTrue(lines.Where(s => !s.StartsWith("#")).All(s => s.Split('\t')[6] != "."));
            File.Delete(filtered);
        }
    }
}