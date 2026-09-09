using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// The VCF entry point annotates a VCF called elsewhere instead of calling variants from reads,
    /// and a bacterial reference requires it because Ensembl Bacteria publishes no known variant
    /// sites for GATK to recalibrate against.
    ///
    /// Two things are worth pinning here. The config key has to be the snake_case name the .smk files
    /// read - the neighbouring "analysisDirectory" key is dead for exactly that reason, and only works
    /// because config.yaml supplies "analysis_directory" itself. And GenerateSpritzCMDArgs has to
    /// re-serialise every option into the in-container invocation: an option omitted there is accepted
    /// on the host and then silently dropped inside the container, which is the classic way a new flag
    /// appears to do nothing.
    /// </summary>
    public class VcfEntryPointTests
    {
        private string _temporaryRoot;

        [SetUp]
        public void SetUp()
        {
            _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-vcf-" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_temporaryRoot);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_temporaryRoot)) Directory.Delete(_temporaryRoot, recursive: true);
        }

        private static SpritzOptions Options() => new()
        {
            Reference = "release-116,homo_sapiens,human,GRCh38",
            SraAccession = "",
            SraAccessionSingleEnd = "",
            Fastq1 = "",
            Fastq1SingleEnd = "",
            AnalyzeVariants = true,
        };

        private string WriteConfig(Action<SpritzOptions> configure)
        {
            SpritzOptions options = Options();
            configure(options);
            var runner = new RunnerEngine(null, _temporaryRoot);
            runner.WriteConfig(options, _temporaryRoot);
            return File.ReadAllText(runner.ConfigFile);
        }

        [Test]
        public void TheVcfIsWrittenUnderTheKeyTheWorkflowReads()
        {
            string config = WriteConfig(o => o.Vcf = "my_variants.vcf");
            Assert.That(config, Does.Contain("vcf: \"my_variants.vcf\""));
        }

        [Test]
        public void AnAbsentVcfIsWrittenEmptySoCheckTreatsItAsUnset()
        {
            // check() in common.smk is len(config[field]) > 0, so "" is what "no VCF" has to look like.
            string config = WriteConfig(o => o.Vcf = null);
            Assert.That(config, Does.Contain("vcf: \"\""));
        }

        [Test]
        public void TheDivisionDefaultsToVertebratesRatherThanBeingAbsent()
        {
            string config = WriteConfig(o => o.Division = null);
            Assert.That(config, Does.Contain("division: \"vertebrates\""));
        }

        [Test]
        public void TheDivisionIsWrittenWhenBacteriaIsAskedFor()
        {
            string config = WriteConfig(o => o.Division = "bacteria");
            Assert.That(config, Does.Contain("division: \"bacteria\""));
        }

        [Test]
        public void TheVcfSurvivesReSerialisationIntoTheContainerInvocation()
        {
            SpritzOptions options = Options();
            options.Vcf = "my_variants.vcf";
            string arguments = SpritzOptionStrings.GenerateSpritzCMDArgs(options);
            Assert.That(arguments, Does.Contain($"--{SpritzOptionStrings.VcfLong}=my_variants.vcf"));
        }

        [Test]
        public void TheDivisionSurvivesReSerialisationIntoTheContainerInvocation()
        {
            SpritzOptions options = Options();
            options.Division = "bacteria";
            string arguments = SpritzOptionStrings.GenerateSpritzCMDArgs(options);
            Assert.That(arguments, Does.Contain($"--{SpritzOptionStrings.DivisionLong}=bacteria"));
        }

        [Test]
        public void AnAbsentVcfAddsNothingToTheContainerInvocation()
        {
            SpritzOptions options = Options();
            options.Vcf = "";
            Assert.That(SpritzOptionStrings.GenerateSpritzCMDArgs(options),
                Does.Not.Contain($"--{SpritzOptionStrings.VcfLong}"));
        }

        [Test]
        public void TheShortVcfOptionDoesNotCollideWithAnotherOption()
        {
            // -v was previously named in stale help text as an analysis flag it has not been for years.
            char[] taken =
            {
                SpritzOptionStrings.AnalysisDirectoryShort, SpritzOptionStrings.AnalyzeVariantsShort,
                SpritzOptionStrings.AnalyzeIsoformsShort, SpritzOptionStrings.QuantifyShort,
                SpritzOptionStrings.AvailableReferencesShort, SpritzOptionStrings.FetchGenomesShort,
                SpritzOptionStrings.AnalysisSetupShort, SpritzOptionStrings.Fastq1Short,
                SpritzOptionStrings.Fastq2Short, SpritzOptionStrings.Fastq1SingleEndShort,
                SpritzOptionStrings.SraAccessionShort, SpritzOptionStrings.SraAccessionSingleEndShort,
                SpritzOptionStrings.ThreadsShort, SpritzOptionStrings.ReferenceShort,
            };
            Assert.That(taken, Does.Not.Contain(SpritzOptionStrings.VcfShort));
            Assert.That(taken, Does.Not.Contain(SpritzOptionStrings.DivisionShort));
            Assert.That(SpritzOptionStrings.VcfShort, Is.Not.EqualTo(SpritzOptionStrings.DivisionShort));
        }

        [Test]
        public void TheStaleHelpTextNoLongerAdvertisesVAsAnAnalysisFlag()
        {
            // It named -v and -w, which have been -b and -c for years; -v now means something else.
            Assert.That(SpritzOptionStrings.InfoSupp, Does.Not.Contain("-v and -w"));
        }

        [TestCase("bacteria", true)]
        [TestCase("Bacteria", true)]
        [TestCase("BACTERIA", true)]
        [TestCase("vertebrates", false)]
        public void BacteriaIsRecognisedWhateverTheCasing(string division, bool expected)
        {
            Assert.That(SpritzOptionStrings.IsBacteria(division), Is.EqualTo(expected));
        }

        [TestCase("vertebrates", true)]
        [TestCase("bacteria", true)]
        [TestCase("Vertebrates", true)]
        [TestCase("fungi", false)]
        [TestCase("", false)]
        public void OnlyTheTwoImplementedDivisionsAreAccepted(string division, bool expected)
        {
            // fungi, protists, plants and metazoa are all real Ensembl Genomes divisions that Spritz
            // does not implement, so they have to be rejected rather than half-work.
            Assert.That(SpritzOptionStrings.IsKnownDivision(division), Is.EqualTo(expected));
        }
    }
}
