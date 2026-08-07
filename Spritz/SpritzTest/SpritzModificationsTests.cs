using NUnit.Framework;
using SpritzModifications;

namespace SpritzTest
{
    public class SpritzModificationsTests
    {
        /// <summary>
        /// Checks that the UniProt ptmlist parses against PSI-MOD into a non-empty modification set, which is
        /// the input to every modification-transfer step in the pipeline.
        ///
        /// Categorised as an external-service test because it downloads. Neither ptmlist.txt nor
        /// PSI-MOD.obo.xml is checked into this repository, and mzLib fetches whatever is missing -
        /// Loaders.LoadPsiMod calls UpdatePsiMod, Loaders.LoadUniprot calls UpdateUniprot - so on a clean
        /// checkout this always goes to the network. The "PSI-MOD database did not exist, writing to disk"
        /// and "Uniprot database did not exist, writing to disk" lines in the test log are that download.
        /// That dependency was previously invisible from the test, which meant an upstream outage looked
        /// like a Spritz regression (compare issues #225 and #233).
        /// </summary>
        [Test]
        [Category("ExternalService")]
        public void GetUniProtMods_ParsesModificationsFromDownloadedPtmList()
        {
            var mods = ProteinAnnotation.GetUniProtMods(TestContext.CurrentContext.TestDirectory);

            Assert.That(mods, Is.Not.Empty,
                "no UniProt modifications were parsed from ptmlist.txt; check the download succeeded and the "
                + "ptmlist format has not changed");
        }
    }
}
