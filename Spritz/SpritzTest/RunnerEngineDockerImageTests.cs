using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Covers how RunnerEngine turns an image name into a command, which is where
    /// RunnerEngine.CurrentVersion is actually consumed: the GUI pulls smithlab/spritz:{CurrentVersion},
    /// so if the compiled version and the published image tag disagree, Spritz pulls an image that does not
    /// exist. That makes this the user-facing end of the version wiring, and it runs offline.
    ///
    /// It also pins the two escape hatches that make local development possible - an explicit tag, and a
    /// non-smithlab image name - which were previously undocumented behaviour of two bare conditionals.
    /// </summary>
    public class RunnerEngineDockerImageTests
    {
        private string _temporaryRoot;
        private RunnerEngine _runner;

        [SetUp]
        public void SetUp()
        {
            // The constructor creates a config directory beside the analysis directory and a resources
            // directory beside that, so it needs somewhere real to write.
            _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-docker-tests-" + Guid.NewGuid().ToString("N"));
            string analysisDirectory = Path.Combine(_temporaryRoot, "analysis");
            Directory.CreateDirectory(analysisDirectory);
            _runner = new RunnerEngine(null, analysisDirectory);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_temporaryRoot))
            {
                Directory.Delete(_temporaryRoot, recursive: true);
            }
        }

        /// <summary>
        /// The default path: the published image is pulled, tagged with the compiled version. This is the
        /// invariant that ties the version to something a user can actually hit.
        /// </summary>
        [Test]
        public void PublishedImageIsTaggedWithTheCompiledVersionAndPulled()
        {
            string command = _runner.GenerateCommandsDry("smithlab/spritz", "SpritzCMD.dll --help");

            Assert.That(command, Does.Contain($"smithlab/spritz:{RunnerEngine.CurrentVersion}"),
                "the published image must be requested at the compiled version, otherwise Spritz pulls a tag "
                + "that release.yml never pushed");
            Assert.That(command, Does.StartWith("podman pull "),
                "the published image should be pulled before it is run, so users get the released container");
        }

        /// <summary>
        /// An explicit tag is honoured verbatim. This is what makes "build a local image and point Spritz at
        /// it" work, and it must not have the compiled version appended to it.
        /// </summary>
        [Test]
        public void ExplicitTagIsUsedVerbatimAndNotPulled()
        {
            string command = _runner.GenerateCommandsDry("spritz:dev", "SpritzCMD.dll --help");

            Assert.That(command, Does.Contain("spritz:dev"));
            Assert.That(command, Does.Not.Contain($"spritz:dev:{RunnerEngine.CurrentVersion}"),
                "an image name that already carries a tag must not have the compiled version appended");
            Assert.That(command, Does.Not.Contain("podman pull"),
                "a locally built image must not be pulled, or the local build would be replaced by whatever "
                + "is in the registry under that name");
        }

        /// <summary>
        /// A local image named without a tag still gets the compiled version appended - it is only the
        /// "smithlab" check that suppresses the pull. Asserted so the two conditions stay distinguishable.
        /// </summary>
        [Test]
        public void UntaggedLocalImageGetsTheCompiledVersionButIsNotPulled()
        {
            string command = _runner.GenerateCommandsDry("spritz", "SpritzCMD.dll --help");

            Assert.That(command, Does.Contain($"spritz:{RunnerEngine.CurrentVersion}"),
                "an untagged name has the compiled version appended, so a local image built as plain "
                + "'spritz' will not be found unless it is tagged with the compiled version");
            Assert.That(command, Does.Not.Contain("podman pull"),
                "only image names containing 'smithlab' are pulled");
        }
    }
}
