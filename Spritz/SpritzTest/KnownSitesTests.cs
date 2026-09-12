using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Which known variant sites base recalibration uses, resolved when the run's config is written.
    ///
    /// It is resolved here rather than in the workflow because deciding it takes a network lookup,
    /// and the workflow re-reads its config on every invocation including dry-runs. The probe is
    /// injectable so these tests never reach Ensembl.
    /// </summary>
    public class KnownSitesTests
    {
        private string _temporaryRoot;

        [SetUp]
        public void SetUp()
        {
            _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-known-" + Guid.NewGuid().ToString("N"));
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
            SraAccession = "SRR629563",
            SraAccessionSingleEnd = "",
            Fastq1 = "",
            Fastq1SingleEnd = "",
            AnalyzeVariants = true,
        };

        /// <summary>Writes a config with the probe stubbed, and returns the file's text.</summary>
        private string WriteConfig(string reference, Func<string, string, bool> probe,
            Action<SpritzOptions> configure = null)
        {
            SpritzOptions options = Options(reference);
            configure?.Invoke(options);
            var runner = new RunnerEngine(null, _temporaryRoot) { EnsemblVariationProbe = probe };
            runner.WriteConfig(options, _temporaryRoot);
            return File.ReadAllText(runner.ConfigFile);
        }

        private static Func<string, string, bool> Answers(bool published) => (_, _) => published;

        private static Func<string, string, bool> Fails =>
            (_, _) => throw new InvalidOperationException("the probe should not have been called");

        [Test]
        public void AutoTakesTheDownloadedSitesWhenEnsemblPublishesThem()
        {
            string config = WriteConfig("release-116,danio_rerio,zebrafish,GRCz11", Answers(true));
            Assert.That(config, Does.Contain("known_sites: \"ensembl\""));
        }

        [Test]
        public void AutoBootstrapsWhenEnsemblPublishesNone()
        {
            // The common case by a wide margin: 340 of 359 species at release 116.
            string config = WriteConfig("release-116,saccharomyces_cerevisiae,yeast,R64-1-1", Answers(false));
            Assert.That(config, Does.Contain("known_sites: \"bootstrap\""));
        }

        [Test]
        public void HumanIsNotProbedBecauseItsSitesComeFromDbSnp()
        {
            // dbSNP is on NCBI, not Ensembl, so asking Ensembl would answer the wrong question.
            string config = WriteConfig("release-116,homo_sapiens,human,GRCh38", Fails);
            Assert.That(config, Does.Contain("known_sites: \"ensembl\""));
        }

        [Test]
        public void BacteriaAreNotProbedBecauseEnsemblBacteriaPublishesNoVariationAtAll()
        {
            string config = WriteConfig(
                "release-63,pseudomonas_aeruginosa_pao1_gca_000006765,pseudomonas,ASM676v1",
                Fails,
                o => o.Division = SpritzOptionStrings.DivisionBacteria);
            Assert.That(config, Does.Contain("known_sites: \"bootstrap\""));
        }

        [Test]
        public void AnExplicitChoiceIsNotProbed()
        {
            string config = WriteConfig("release-116,danio_rerio,zebrafish,GRCz11", Fails,
                o => o.KnownSites = EnsemblVariation.Bootstrap);
            Assert.That(config, Does.Contain("known_sites: \"bootstrap\""));
        }

        [Test]
        public void AnExplicitEnsemblChoiceSurvivesAProbeThatWouldSayOtherwise()
        {
            string config = WriteConfig("release-116,saccharomyces_cerevisiae,yeast,R64-1-1",
                Answers(false), o => o.KnownSites = EnsemblVariation.Ensembl);
            Assert.That(config, Does.Contain("known_sites: \"ensembl\""));
        }

        [Test]
        public void TheWorkflowOnlyEverSeesAResolvedMode()
        {
            // common.smk accepts "ensembl" or "bootstrap" and rejects anything else, so "auto" must
            // never reach it.
            string config = WriteConfig("release-116,danio_rerio,zebrafish,GRCz11", Answers(true),
                o => o.KnownSites = EnsemblVariation.Auto);
            // Asserted positively: a bare Does.Not.Contain passed even with the whole known_sites
            // block deleted, which is the opposite of pinning an invariant.
            Assert.That(config, Does.Contain("known_sites: \"ensembl\""));
            Assert.That(config, Does.Not.Contain("known_sites: \"auto\""));
        }

        [Test]
        public void TheModeSurvivesReSerialisationIntoTheContainerInvocation()
        {
            SpritzOptions options = Options("release-116,danio_rerio,zebrafish,GRCz11");
            options.KnownSites = EnsemblVariation.Bootstrap;
            Assert.That(SpritzOptionStrings.GenerateSpritzCMDArgs(options),
                Does.Contain($"--{SpritzOptionStrings.KnownSitesLong}=bootstrap"));
        }

        [TestCase("auto", true)]
        [TestCase("ensembl", true)]
        [TestCase("bootstrap", true)]
        [TestCase("Bootstrap", true)]
        [TestCase("vqsr", false)]
        [TestCase("", false)]
        public void OnlyTheThreeModesAreAccepted(string mode, bool expected)
        {
            Assert.That(EnsemblVariation.IsKnownMode(mode), Is.EqualTo(expected));
        }

        [Test]
        public void TheProbeAsksForTheFilenameTheDownloadRuleActuallyFetches()
        {
            // Only the lowercase form resolves for any species today, so it is checked first;
            // the rule itself tries the capitalised name first and falls back. Checking the
            // directory alone, or only the capitalised name, would answer the wrong question.
            string[] urls = EnsemblVariation.CandidateUrls("116", "danio_rerio");
            Assert.That(urls[0], Is.EqualTo(
                "https://ftp.ensembl.org/pub/release-116/variation/vcf/danio_rerio/danio_rerio.vcf.gz"));
            Assert.That(urls[1], Is.EqualTo(
                "https://ftp.ensembl.org/pub/release-116/variation/vcf/danio_rerio/Danio_rerio.vcf.gz"));
        }

        [Test]
        public void TheShortOptionDoesNotCollide()
        {
            char[] taken =
            {
                SpritzOptionStrings.AnalysisDirectoryShort, SpritzOptionStrings.AnalyzeVariantsShort,
                SpritzOptionStrings.AnalyzeIsoformsShort, SpritzOptionStrings.QuantifyShort,
                SpritzOptionStrings.AvailableReferencesShort, SpritzOptionStrings.FetchGenomesShort,
                SpritzOptionStrings.AnalysisSetupShort, SpritzOptionStrings.Fastq1Short,
                SpritzOptionStrings.Fastq2Short, SpritzOptionStrings.Fastq1SingleEndShort,
                SpritzOptionStrings.SraAccessionShort, SpritzOptionStrings.SraAccessionSingleEndShort,
                SpritzOptionStrings.ThreadsShort, SpritzOptionStrings.ReferenceShort,
                SpritzOptionStrings.VcfShort, SpritzOptionStrings.DivisionShort,
            };
            Assert.That(taken, Does.Not.Contain(SpritzOptionStrings.KnownSitesShort));
        }
    }
}
