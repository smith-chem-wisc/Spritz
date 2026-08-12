using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Covers the three runtimes Spritz can drive. Podman is the default because it is Apache-2.0 and
    /// installs without Docker Desktop; Docker stays for users who already have it; Apptainer exists
    /// for HPC, where a root daemon is usually not permitted.
    /// </summary>
    public class ContainerRuntimeTests
    {
        private string _temporaryRoot;
        private RunnerEngine _runner;

        [SetUp]
        public void SetUp()
        {
            _temporaryRoot = Path.Combine(Path.GetTempPath(), "spritz-runtime-tests-" + Guid.NewGuid().ToString("N"));
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

        private string Command(ContainerRuntime runtime)
        {
            _runner.Runtime = runtime;
            return _runner.GenerateCommandsDry("smithlab/spritz", "SpritzCMD.dll --help");
        }

        [Test]
        public void PodmanIsTheDefault()
        {
            Assert.That(new RunnerEngine(null, Path.GetTempPath()).Runtime, Is.EqualTo(ContainerRuntime.Podman));
        }

        [TestCase(ContainerRuntime.Podman, "podman")]
        [TestCase(ContainerRuntime.Docker, "docker")]
        public void DockerCompatibleRuntimesDifferOnlyByExecutable(ContainerRuntime runtime, string exe)
        {
            string command = Command(runtime);
            Assert.That(command, Does.StartWith($"{exe} pull "));
            Assert.That(command, Does.Contain($"{exe} run --rm -i -t --user=root --name "));
            Assert.That(command, Does.Contain($"{exe} stop "));
        }

        [Test]
        public void PodmanAndDockerProduceTheSameCommandApartFromTheExecutable()
        {
            Assert.That(Command(ContainerRuntime.Podman).Replace("podman", "docker"),
                Is.EqualTo(Command(ContainerRuntime.Docker)),
                "these two share a builder, so any divergence is a bug rather than a decision");
        }

        /// <summary>
        /// Apptainer has no daemon, so naming or stopping a container is meaningless, and it binds
        /// rather than mounts. It reads the same published image through a docker:// URI.
        /// </summary>
        [Test]
        public void ApptainerBindsAndUsesADockerUri()
        {
            string command = Command(ContainerRuntime.Apptainer);

            Assert.That(command, Does.StartWith("apptainer run "));
            Assert.That(command, Does.Contain("docker://smithlab/spritz:"));
            Assert.That(command, Does.Contain("--bind "));
            Assert.That(command, Does.Not.Contain("-v "), "apptainer binds, it does not mount");
            Assert.That(command, Does.Not.Contain("--name"), "there is no daemon to name a container in");
            Assert.That(command, Does.Not.Contain("stop"), "and nothing to stop afterwards");
        }

        /// <summary>
        /// The image is read-only under Apptainer and the workflow writes inside it, so without this
        /// the run fails partway through rather than at startup.
        /// </summary>
        [Test]
        public void ApptainerMakesTheImageWritable()
        {
            Assert.That(Command(ContainerRuntime.Apptainer), Does.Contain("--writable-tmpfs"));
        }

        [Test]
        public void EveryRuntimeMountsBothWorkingDirectories()
        {
            foreach (ContainerRuntime runtime in Enum.GetValues<ContainerRuntime>())
            {
                string command = Command(runtime);
                Assert.That(command, Does.Contain("/app/spritz/results/"), runtime.ToString());
                Assert.That(command, Does.Contain("/app/spritz/resources"), runtime.ToString());
            }
        }

        [Test]
        public void TopIsEmptyForApptainerBecauseThereIsNoDaemon()
        {
            _runner.Runtime = ContainerRuntime.Apptainer;
            Assert.That(_runner.GenerateTopComand(), Is.Empty);

            _runner.Runtime = ContainerRuntime.Podman;
            Assert.That(_runner.GenerateTopComand(), Does.StartWith("podman container top "));
        }

        [TestCase("podman", ContainerRuntime.Podman)]
        [TestCase("Docker", ContainerRuntime.Docker)]
        [TestCase("apptainer", ContainerRuntime.Apptainer)]
        [TestCase("singularity", ContainerRuntime.Apptainer)]
        [TestCase("", ContainerRuntime.Podman)]
        [TestCase(null, ContainerRuntime.Podman)]
        [TestCase("nonsense", ContainerRuntime.Podman)]
        public void RuntimeNamesParseLenientlyAndDefaultToPodman(string name, ContainerRuntime expected)
        {
            Assert.That(ContainerRuntimes.Parse(name), Is.EqualTo(expected));
        }
    }
}
