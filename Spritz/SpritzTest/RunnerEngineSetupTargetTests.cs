using NUnit.Framework;
using SpritzBackend;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Covers the snakemake target used by "Analysis Setup". The rule is named setup and writes
    /// ../resources/setup.txt, so the previous " setup.txt" matched neither and snakemake refused
    /// with MissingRuleException before any work started.
    /// </summary>
    public class RunnerEngineSetupTargetTests
    {
        private static string Command(bool setup) =>
            new RunnerEngine(null, Path.GetTempPath())
                .GenerateSnakemakeCommand(new SpritzOptions { Threads = 4 }, setup);

        [Test]
        public void SetupRequestsTheRuleByName()
        {
            Assert.That(Command(setup: true), Does.EndWith(" setup"));
        }

        [Test]
        public void SetupDoesNotAskForAFilenameThatMatchesNoRule()
        {
            Assert.That(Command(setup: true), Does.Not.Contain("setup.txt"));
        }

        [Test]
        public void ANormalRunNamesNoTargetSoTheDefaultRuleIsUsed()
        {
            Assert.That(Command(setup: false), Does.Not.Contain("setup"));
        }

        [Test]
        public void BothFormsStillCarryTheCoreFlags()
        {
            foreach (bool setup in new[] { true, false })
            {
                string command = Command(setup);
                Assert.That(command, Does.StartWith("snakemake "), command);
                Assert.That(command, Does.Contain("--use-conda"), command);
                Assert.That(command, Does.Contain("--configfile"), command);
            }
        }
    }
}
