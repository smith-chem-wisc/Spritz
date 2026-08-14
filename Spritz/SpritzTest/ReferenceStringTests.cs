using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// The reference string is a whole line copied out of genomes.csv, e.g.
    /// "release-116,homo_sapiens,human,GRCh38". WriteConfig splits it into four fields and writes them
    /// into config.yaml, and everything downstream keys off those values, so a wrong one is worth
    /// reporting precisely rather than as an index error or a type name.
    /// </summary>
    public class ReferenceStringTests
    {
        private string _temporaryRoot;

        [SetUp]
        public void SetUp()
        {
            _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-reference-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_temporaryRoot);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_temporaryRoot)) Directory.Delete(_temporaryRoot, recursive: true);
        }

        private static SpritzOptions Options(string reference) => new()
        {
            Reference = reference,
            SraAccession = "",
            SraAccessionSingleEnd = "",
            Fastq1 = "",
            Fastq1SingleEnd = "",
            AnalyzeVariants = true,
        };

        private string WriteConfig(string reference)
        {
            var runner = new RunnerEngine(null, _temporaryRoot);
            runner.WriteConfig(Options(reference), _temporaryRoot);
            return File.ReadAllText(runner.ConfigFile);
        }

        [Test]
        public void AWholeGenomesCsvLineIsWrittenToTheConfig()
        {
            string config = WriteConfig("release-116,homo_sapiens,human,GRCh38");

            Assert.That(config, Contains.Substring("release: \"116\""), "the release- prefix is stripped");
            Assert.That(config, Contains.Substring("species: \"Homo_sapiens\""), "species is capitalised");
            Assert.That(config, Contains.Substring("organism: \"human\""));
            Assert.That(config, Contains.Substring("genome: \"GRCh38\""));
        }

        /// <summary>
        /// The message used to interpolate the string[] that Split returned, so it read
        /// 'the reference string "System.String[]" does not have four comma-separated elements' -
        /// naming a type instead of what the user typed.
        /// </summary>
        [Test]
        [TestCase("release-116,homo_sapiens,human", TestName = "TooFewFields")]
        [TestCase("release-116,homo_sapiens,human,GRCh38,extra", TestName = "TooManyFields")]
        [TestCase("GRCh38", TestName = "OneField")]
        public void AReferenceWithoutFourFieldsNamesWhatWasTyped(string reference)
        {
            var thrown = Assert.Throws<SpritzException>(() => WriteConfig(reference));

            Assert.That(thrown.Message, Contains.Substring(reference));
            Assert.That(thrown.Message, Does.Not.Contain("System.String[]"));
            Assert.That(thrown.Message, Contains.Substring("genomes.csv"), "say where a correct one comes from");
        }

        /// <summary>
        /// An empty field used to reach further in and fail as an index error, or produce an empty value
        /// in config.yaml that only failed much later in the workflow.
        /// </summary>
        [Test]
        [TestCase(",homo_sapiens,human,GRCh38", 1)]
        [TestCase("release-116,,human,GRCh38", 2)]
        [TestCase("release-116,homo_sapiens,,GRCh38", 3)]
        [TestCase("release-116,homo_sapiens,human,", 4)]
        public void AnEmptyFieldIsReportedByPosition(string reference, int position)
        {
            var thrown = Assert.Throws<SpritzException>(() => WriteConfig(reference));

            Assert.That(thrown.Message, Contains.Substring($"element {position}"));
            Assert.That(thrown.Message, Contains.Substring(reference));
        }

        /// <summary>
        /// The release number was read as releaseStr[8..], a bare index that happened to be the length of
        /// "release-". Anything shorter threw ArgumentOutOfRangeException, which says nothing about which
        /// field was wrong or what it should look like.
        /// </summary>
        [Test]
        [TestCase("release-116,homo_sapiens,human,GRCh38", "116", TestName = "PrefixedRelease")]
        [TestCase("116,homo_sapiens,human,GRCh38", "116", TestName = "BareNumber")]
        [TestCase("RELEASE-116,homo_sapiens,human,GRCh38", "116", TestName = "PrefixCasingIgnored")]
        [TestCase("release-81,bos_taurus,bovine,UMD3.1", "81", TestName = "OldestLineInGenomesCsv")]
        public void TheReleaseNumberIsReadFromEitherForm(string reference, string expected)
        {
            Assert.That(WriteConfig(reference), Contains.Substring($"release: \"{expected}\""));
        }

        [Test]
        [TestCase("rel,homo_sapiens,human,GRCh38", TestName = "ShorterThanThePrefix")]
        [TestCase("release-,homo_sapiens,human,GRCh38", TestName = "PrefixWithNoNumber")]
        [TestCase("release-abc,homo_sapiens,human,GRCh38", TestName = "NotANumber")]
        [TestCase("release-11.6,homo_sapiens,human,GRCh38", TestName = "NotAnInteger")]
        public void SomethingThatIsNotAReleaseIsRejectedByName(string reference)
        {
            var thrown = Assert.Throws<SpritzException>(() => WriteConfig(reference));

            Assert.That(thrown.Message, Contains.Substring("release"),
                "name the field, not just the position");
            Assert.That(thrown.Message, Contains.Substring("release-116"), "show a correct example");
        }

        /// <summary>
        /// Capitalisation used to go through string.ToUpper(), which is culture-sensitive: under a Turkish
        /// locale a leading "i" becomes a dotted capital and the species no longer matches the Ensembl
        /// directory name.
        /// </summary>
        [Test]
        public void SpeciesCapitalisationDoesNotDependOnTheCurrentCulture()
        {
            var previous = System.Globalization.CultureInfo.CurrentCulture;
            try
            {
                System.Globalization.CultureInfo.CurrentCulture = new System.Globalization.CultureInfo("tr-TR");
                Assert.That(WriteConfig("release-116,ictalurus_punctatus,channel catfish,IpCoco_1.2"),
                    Contains.Substring("species: \"Ictalurus_punctatus\""));
            }
            finally
            {
                System.Globalization.CultureInfo.CurrentCulture = previous;
            }
        }

        /// <summary>
        /// A single-character species field would have thrown on the [1..] slice.
        /// </summary>
        [Test]
        public void ASingleCharacterSpeciesIsNotAnIndexError()
        {
            Assert.That(WriteConfig("release-116,x,thing,ASM1"), Contains.Substring("species: \"X\""));
        }
    }
}
