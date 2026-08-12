using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;
using System.Linq;

namespace SpritzTest
{
    /// <summary>
    /// Covers the three fragile things in how the container command used to be built: a container name
    /// from String.GetHashCode, a pull decision from a substring search, and paths escaped into a single
    /// shell string.
    /// </summary>
    public class RunnerEngineHardeningTests
    {
        private string _temporaryRoot;

        private RunnerEngine Runner(string analysisSubdirectory = "analysis")
        {
            string directory = Path.Combine(_temporaryRoot, analysisSubdirectory);
            Directory.CreateDirectory(directory);
            return new RunnerEngine(null, directory);
        }

        [SetUp]
        public void SetUp() => _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-hardening-" + Guid.NewGuid().ToString("N"));

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_temporaryRoot)) Directory.Delete(_temporaryRoot, recursive: true);
        }

        // ---- container name ------------------------------------------------------

        /// <summary>
        /// String.GetHashCode is randomised per process in .NET Core, so a second process could not
        /// work out the name of a container this one started - which is what `docker stop` needs.
        /// </summary>
        [Test]
        public void ContainerNameIsReadableAndCarriesNoHash()
        {
            string name = Runner().SpritzContainerName;

            Assert.That(name, Does.StartWith("spritz-"));
            Assert.That(name, Does.Match(@"^spritz-\d{4}-\d{2}-\d{2}-\d{2}-\d{2}-\d{2}$"),
                $"expected a timestamp rather than a hash, got '{name}'");
        }

        [Test]
        public void ContainerNameIsLegalForDockerAndPodman()
        {
            // Both require [a-zA-Z0-9][a-zA-Z0-9_.-]* - a negative hash would have started with '-'.
            Assert.That(Runner().SpritzContainerName, Does.Match(@"^[a-zA-Z0-9][a-zA-Z0-9_.-]*$"));
        }

        [Test]
        public void TheStopCommandNamesTheContainerThatWasStarted()
        {
            var runner = Runner();
            string command = runner.GenerateCommandsDry("smithlab/spritz", "SpritzCMD.dll --help");

            Assert.That(command.Split(runner.SpritzContainerName).Length - 1, Is.EqualTo(2),
                "the name should appear twice - once to create, once to stop - and be the same both times");
        }

        // ---- pull decision -------------------------------------------------------

        [TestCase("smithlab/spritz", true)]
        [TestCase("docker.io/smithlab/spritz", true)]
        [TestCase("ghcr.io/smith-chem-wisc/spritz", true)]
        [TestCase("SmithLab/spritz", true)]
        [TestCase("spritz", false)]
        [TestCase("spritz:dev", false)]
        [TestCase("myfork/spritz", false)]
        public void OnlyImagesFromAPublishedRegistryArePulled(string image, bool expected)
        {
            Assert.That(Runner().ShouldPull(image), Is.EqualTo(expected));
        }

        /// <summary>The substring check pulled this and overwrote a local build with the published one.</summary>
        [TestCase("smithlab-test")]
        [TestCase("my-smithlab-experiment")]
        [TestCase("localhost:5000/smithlab-copy")]
        public void ALocalImageThatMerelyMentionsTheRegistryIsNotPulled(string image)
        {
            Assert.That(Runner().ShouldPull(image), Is.False, image);
        }

        [Test]
        public void AForkCanOptIntoPullingItsOwnImage()
        {
            var runner = Runner();
            Assert.That(runner.ShouldPull("myfork/spritz"), Is.False);
            runner.AlwaysPull = true;
            Assert.That(runner.ShouldPull("myfork/spritz"), Is.True);
        }

        // ---- arguments rather than a shell string --------------------------------

        [Test]
        public void ArgumentsAreNotQuotedOrEscaped()
        {
            var (executable, arguments) = Runner().GenerateRunArguments("smithlab/spritz", "SpritzCMD.dll --help");

            Assert.That(executable, Is.EqualTo("podman"));
            Assert.That(arguments.Any(a => a.Contains('"')), Is.False, "nothing should need quoting");
            Assert.That(arguments, Does.Contain("--name"));
        }

        /// <summary>
        /// The whole point of the argument list: a directory containing a space arrives as one argument
        /// instead of being split by a shell.
        /// </summary>
        [Test]
        public void APathWithASpaceStaysOneArgument()
        {
            var (_, arguments) = Runner("analysis with spaces").GenerateRunArguments("smithlab/spritz", "SpritzCMD.dll --help");

            string mount = arguments.SingleOrDefault(a => a.Contains("analysis with spaces"));
            Assert.That(mount, Is.Not.Null, string.Join(" | ", arguments));
            Assert.That(mount, Does.EndWith(":/app/spritz/results/"));
        }

        [Test]
        public void ApptainerBindsRatherThanMountsInTheArgumentListToo()
        {
            var runner = Runner();
            runner.Runtime = ContainerRuntime.Apptainer;
            var (executable, arguments) = runner.GenerateRunArguments("smithlab/spritz", "SpritzCMD.dll --help");

            Assert.That(executable, Is.EqualTo("apptainer"));
            Assert.That(arguments, Does.Contain("--bind"));
            Assert.That(arguments, Does.Not.Contain("-v"));
            Assert.That(arguments.Any(a => a.StartsWith("docker://")), Is.True);
        }
    }
}
