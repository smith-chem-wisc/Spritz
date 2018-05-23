using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class TestSequenceSimilarityMethods
    {
        [Test]
        public void test_modification_transfer_exact_sequence_match()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "xml2.xml"), true, DecoyType.None, null, false, null, out Dictionary<string, Modification> un);
            List<Protein> destination = new List<Protein> {
                new Protein("MKTCYYELLGVETHASDLELKKAYRKKALQYHPDKNPDNVEEATQKFAVIRAAYEVLSDPQERAWYDSHKEQILNDTPPSTDDYYDYEVDATVTGVTTDELLLFFNSALYTKIDNSAAGIYQIAGKIFAKLAKDEILSGKRLGKFSEYQDDVFEQDINSIGYLKACDNFINKTDKLLYPLFGYSPTDYEYLKHFYKTWSAFNTLKSFSWKDEYMYSKNYDRRTKREVNRRNEKARQQARNEYNKTVKRFVVFIKKLDKRMKEGAKIAEEQRKLKEQQRKNELNNRRKFGNDNNDEEKFHLQSWQTVKEENWDELEKVYDNFGEFENSKNDKEGEVLIYECFICNKTFKSEKQLKNHINTKLHKKNMEEIRKEMEEENITLGLDNLSDLEKFDSADESVKEKEDIDLQALQAELAEIERKLAESSSEDESEDDNLNIEMDIEVEDVSSDENVHVNTKNKKKRKKKKKAKVDTETEESESFDDTKDKRSNELDDLLASLGDKGLQTDDDEDWSTKAKKKKGKQPKKNSKSTKSTPSLSTLPSSMSPTSAIEVCTTCGESFDSRNKLFNHVKIAGHAAVKNVVKRKKVKTKRI",
                    "") };

            Assert.AreEqual(ok[0].BaseSequence, destination[0].BaseSequence);
            List<Protein> new_proteins = ProteinAnnotation.TransferModifications(ok, destination);

            Assert.AreEqual(1, new_proteins.Count);
            Assert.AreEqual(ok[0].BaseSequence, new_proteins[0].BaseSequence);
            Assert.AreEqual(ok[0].OneBasedPossibleLocalizedModifications, new_proteins[0].OneBasedPossibleLocalizedModifications);
            Assert.IsTrue(new_proteins[0].OneBasedPossibleLocalizedModifications.Keys.Count == 2);
            Assert.IsTrue(ok[0].DatabaseReferences.All(x => new_proteins[0].DatabaseReferences.Contains(x)));
            Assert.IsTrue(ok[0].ProteolysisProducts.All(x => new_proteins[0].ProteolysisProducts.Contains(x)));
        }
    }
}
