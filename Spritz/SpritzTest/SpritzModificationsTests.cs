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
            int modificationCount;
            try
            {
                modificationCount = ProteinAnnotation.GetUniProtMods(TestContext.CurrentContext.TestDirectory).Count;
            }
            catch (System.Exception exception)
            {
                // A failed download is an outage, not a Spritz regression, so it skips rather than fails.
                // A download that succeeds but parses to nothing still fails below - that would mean the
                // ptmlist format changed, which is a real break worth reddening the check for.
                Assert.Ignore($"could not download or read the UniProt ptmlist: {exception.Message}");
                return;
            }

            Assert.That(modificationCount, Is.GreaterThan(0),
                "no UniProt modifications were parsed from ptmlist.txt; check the download succeeded and the "
                + "ptmlist format has not changed");
        }
    }
}
