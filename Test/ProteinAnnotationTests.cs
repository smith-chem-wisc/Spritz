using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using WorkflowLayer;

namespace Test
{
    [TestFixture]
    public class ProteinAnnotationTests
    {
        /// <summary>
        /// This really should be in Proteogenomics, but it's failing with weird invalid program detected errors.
        /// </summary>
        [Test]
        public void ProteinAnnTransferExactSequenceMatchMods()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "xml2.xml"), true, DecoyType.None, null, false, null, out Dictionary<string, Modification> un);
            List<Protein> destination = new List<Protein> {
                new Protein("MKTCYYELLGVETHASDLELKKAYRKKALQYHPDKNPDNVEEATQKFAVIRAAYEVLSDPQERAWYDSHKEQILNDTPPSTDDYYDYEVDATVTGVTTDELLLFFNSALYTKIDNSAAGIYQIAGKIFAKLAKDEILSGKRLGKFSEYQDDVFEQDINSIGYLKACDNFINKTDKLLYPLFGYSPTDYEYLKHFYKTWSAFNTLKSFSWKDEYMYSKNYDRRTKREVNRRNEKARQQARNEYNKTVKRFVVFIKKLDKRMKEGAKIAEEQRKLKEQQRKNELNNRRKFGNDNNDEEKFHLQSWQTVKEENWDELEKVYDNFGEFENSKNDKEGEVLIYECFICNKTFKSEKQLKNHINTKLHKKNMEEIRKEMEEENITLGLDNLSDLEKFDSADESVKEKEDIDLQALQAELAEIERKLAESSSEDESEDDNLNIEMDIEVEDVSSDENVHVNTKNKKKRKKKKKAKVDTETEESESFDDTKDKRSNELDDLLASLGDKGLQTDDDEDWSTKAKKKKGKQPKKNSKSTKSTPSLSTLPSSMSPTSAIEVCTTCGESFDSRNKLFNHVKIAGHAAVKNVVKRKKVKTKRI",
                    "") };

            Assert.AreEqual(ok[0].BaseSequence, destination[0].BaseSequence);
            List<Protein> newProteins = ProteinAnnotation.CombineAndAnnotateProteins(ok, destination);

            Assert.AreEqual(1, newProteins.Count);
            Assert.AreEqual(ok[0].BaseSequence, newProteins[0].BaseSequence);
            Assert.AreEqual(ok[0].OneBasedPossibleLocalizedModifications, newProteins[0].OneBasedPossibleLocalizedModifications);
            Assert.IsTrue(newProteins[0].OneBasedPossibleLocalizedModifications.Keys.Count == 2);
            Assert.IsTrue(ok[0].DatabaseReferences.All(x => newProteins[0].DatabaseReferences.Contains(x)));
            Assert.IsTrue(ok[0].ProteolysisProducts.All(x => newProteins[0].ProteolysisProducts.Contains(x)));
        }

        [Test]
        public void ProteinAnnCombineSequenceEntries()
        {
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "xml2.xml"), true, DecoyType.None, null, false, null, out Dictionary<string, Modification> un);
            List<Protein> destination = new List<Protein> {
                new Protein("MKTCYYELLGVETHASDLELKKAYRKKALQYHPDKNPDNVEEATQKFAVIRAAYEVLSDPQERAWYDSHKEQILNDTPPSTDDYYDYEVDATVTGVTTDELLLFFNSALYTKIDNSAAGIYQIAGKIFAKLAKDEILSGKRLGKFSEYQDDVFEQDINSIGYLKACDNFINKTDKLLYPLFGYSPTDYEYLKHFYKTWSAFNTLKSFSWKDEYMYSKNYDRRTKREVNRRNEKARQQARNEYNKTVKRFVVFIKKLDKRMKEGAKIAEEQRKLKEQQRKNELNNRRKFGNDNNDEEKFHLQSWQTVKEENWDELEKVYDNFGEFENSKNDKEGEVLIYECFICNKTFKSEKQLKNHINTKLHKKNMEEIRKEMEEENITLGLDNLSDLEKFDSADESVKEKEDIDLQALQAELAEIERKLAESSSEDESEDDNLNIEMDIEVEDVSSDENVHVNTKNKKKRKKKKKAKVDTETEESESFDDTKDKRSNELDDLLASLGDKGLQTDDDEDWSTKAKKKKGKQPKKNSKSTKSTPSLSTLPSSMSPTSAIEVCTTCGESFDSRNKLFNHVKIAGHAAVKNVVKRKKVKTKRI",
                    "Acc1", organism: "Homo sapiens", gene_names: new List<Tuple<string, string>>{ new Tuple<string, string>( "primary", "gene1" ) },
                    name: "name1", full_name: "fullname1"),
                new Protein("MKTCYYELLGVETHASDLELKKAYRKKALQYHPDKNPDNVEEATQKFAVIRAAYEVLSDPQERAWYDSHKEQILNDTPPSTDDYYDYEVDATVTGVTTDELLLFFNSALYTKIDNSAAGIYQIAGKIFAKLAKDEILSGKRLGKFSEYQDDVFEQDINSIGYLKACDNFINKTDKLLYPLFGYSPTDYEYLKHFYKTWSAFNTLKSFSWKDEYMYSKNYDRRTKREVNRRNEKARQQARNEYNKTVKRFVVFIKKLDKRMKEGAKIAEEQRKLKEQQRKNELNNRRKFGNDNNDEEKFHLQSWQTVKEENWDELEKVYDNFGEFENSKNDKEGEVLIYECFICNKTFKSEKQLKNHINTKLHKKNMEEIRKEMEEENITLGLDNLSDLEKFDSADESVKEKEDIDLQALQAELAEIERKLAESSSEDESEDDNLNIEMDIEVEDVSSDENVHVNTKNKKKRKKKKKAKVDTETEESESFDDTKDKRSNELDDLLASLGDKGLQTDDDEDWSTKAKKKKGKQPKKNSKSTKSTPSLSTLPSSMSPTSAIEVCTTCGESFDSRNKLFNHVKIAGHAAVKNVVKRKKVKTKRI",
                    "Acc2", organism: "Homo sapiens", gene_names: new List<Tuple<string, string>>{ new Tuple<string, string>( "primary", "gene2" ) },
                    name: "name2", full_name: "fullname2"),
                new Protein("MNOTTHESAMESEQ",
                    "Acc2", organism: "Homo sapiens", gene_names: new List<Tuple<string, string>>{ new Tuple<string, string>( "primary", "gene2" ) },
                    name: "name2", full_name: "fullname2"),
            };

            List<Protein> newProteins = ProteinAnnotation.CombineAndAnnotateProteins(ok, destination);

            Assert.AreEqual(2, newProteins.Count); // two were combined
            Assert.IsTrue(newProteins.Any(p => p.Name.Contains(destination[0].Name) && p.Name.Contains(destination[1].Name)));
        }
    }
}