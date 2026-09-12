using NUnit.Framework;
using SpritzBackend;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Covers the --resources limit that serializes the rules running SpritzModifications on the
    /// UniProt xml. A rule's `resources: uniprot_temp=1` constrains scheduling only when a limit for
    /// that name is given on the command line; without the flag snakemake treats it as unlimited and
    /// the four rules run concurrently again, colliding on the fixed temp.xml mzLib decompresses to.
    /// So the token in proteogenomics.smk and the limit here are one fix in two halves, and dropping
    /// either half restores the crash silently.
    /// </summary>
    public class RunnerEngineResourceLimitTests
    {
        private static string Command(bool setup) =>
            new RunnerEngine(null, Path.GetTempPath())
                .GenerateSnakemakeCommand(new SpritzOptions { Threads = 4 }, setup);

        [Test]
        public void BothFormsLimitTheUniProtTempResource()
        {
            foreach (bool setup in new[] { true, false })
            {
                string command = Command(setup);
                Assert.That(command, Does.Contain("--resources uniprot_temp=1"), command);
            }
        }

        [Test]
        public void TheLimitIsOneSoTheRulesCannotOverlap()
        {
            // A limit above 1 would let two of them run together, which is the whole failure.
            Assert.That(Command(setup: false), Does.Not.Match(@"uniprot_temp=(?!1\b)\d+"));
        }

        [Test]
        public void ThreadsStillDriveJobParallelismRatherThanTheTokenReplacingIt()
        {
            // The token must not have been mistaken for a global serialization: everything else
            // still runs -j wide.
            Assert.That(Command(setup: false), Does.Contain("-j 4"));
        }
    }
}
