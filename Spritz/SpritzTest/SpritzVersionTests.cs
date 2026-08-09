using NUnit.Framework;
using SpritzBackend;
using System.IO;
using System.Linq;

namespace SpritzTest
{
    public class SpritzVersionTests
    {
        /// <summary>
        /// Guards against publishing a release newer than the version compiled into the code, i.e. against
        /// forgetting to bump RunnerEngine.CurrentVersion.
        ///
        /// Queries the GitHub releases API, so it is categorised as an external-service test. It cannot pass
        /// offline, and the unauthenticated GitHub API is rate limited per IP - which CI runners share - so it
        /// is flaky by nature. It is excluded from the required CI run and executed in a separate
        /// non-blocking job instead.
        ///
        /// When the API cannot answer, this reports as SKIPPED rather than failed. The point of the test is to
        /// nag about bumping RunnerEngine.CurrentVersion, and that nag only carries if a red result means the
        /// versions really disagree. Rate-limit noise reddening the check trains everyone to ignore it.
        /// </summary>
        [Test]
        [Category("ExternalService")]
        public void PublishedReleaseIsNotNewerThanCompiledVersion()
        {
            SpritzVersion version = new();
            try
            {
                version.GetVersionNumbersFromWeb();
            }
            catch (System.Exception exception)
            {
                Assert.Ignore($"could not reach the GitHub releases API: {exception.Message}");
            }

            // Distinguish "the API told us nothing" from "the versions disagree", and skip rather than fail in
            // the first case. A rate-limit reply is valid JSON with no tag_name, which leaves this null.
            if (string.IsNullOrEmpty(version.NewestKnownVersion))
            {
                Assert.Ignore("the GitHub releases API returned no tag_name, most likely a rate limit or an "
                    + "outage; skipped rather than failed because this is not a version mismatch");
            }

            // SpritzVersion.IsVersionLower(null) throws NullReferenceException rather than returning a
            // value: GetVersionNumber catches FormatException but not a null input. The no-installer case is
            // therefore handled here rather than being passed through it.
            if (version.NewestKnownVersionWithMsi is null)
            {
                Assert.That(version.NewestKnownVersion, Is.EqualTo(RunnerEngine.CurrentVersion),
                    $"newest published release is {version.NewestKnownVersion} and none carries an installer, "
                    + $"but the code reports {RunnerEngine.CurrentVersion}");
                return;
            }

            Assert.That(SpritzVersion.IsVersionLower(version.NewestKnownVersionWithMsi), Is.True,
                $"newest published release with an installer is {version.NewestKnownVersionWithMsi}, which is newer "
                + $"than the compiled version {RunnerEngine.CurrentVersion}; bump RunnerEngine.CurrentVersion");
        }

        /// <summary>
        /// The installer's Product/@Version and RunnerEngine.CurrentVersion must agree, otherwise a release
        /// ships an installer labelled with a different version than the application reports. Runs offline.
        /// </summary>
        [Test]
        public void InstallerVersionMatchesCompiledVersion()
        {
            FileInfo productWxs = FindProductWxs();

            string productLine = File.ReadLines(productWxs.FullName)
                .FirstOrDefault(line => line.TrimStart().StartsWith("<Product"));
            Assert.That(productLine, Is.Not.Null,
                $"no <Product ...> element found in {productWxs.FullName}");

            string versionAttribute = productLine.Split(' ')
                .FirstOrDefault(token => token.StartsWith("Version"));
            Assert.That(versionAttribute, Is.Not.Null,
                $"the <Product ...> element in {productWxs.FullName} has no Version attribute: {productLine.Trim()}");

            string installerVersion = versionAttribute.Split('=')[1].Trim('"');
            Assert.That(installerVersion, Is.EqualTo(RunnerEngine.CurrentVersion),
                $"{productWxs.Name} declares Version={installerVersion} but RunnerEngine.CurrentVersion is "
                + $"{RunnerEngine.CurrentVersion}; these must be bumped together");
        }

        /// <summary>
        /// Locates SpritzInstaller/Product.wxs by walking up from the test assembly.
        ///
        /// This was previously a fixed "../../../../../" hop from Environment.CurrentDirectory, which encoded
        /// the depth of the build output directory - bin/{Platform}/{Configuration}/{TargetFramework} - and so
        /// would silently start looking in the wrong place as soon as the platform or framework segment
        /// changed. Searching upward for a known marker does not care how deep the output happens to be.
        /// </summary>
        private static FileInfo FindProductWxs()
        {
            string startedAt = TestContext.CurrentContext.TestDirectory;
            for (DirectoryInfo directory = new(startedAt); directory is not null; directory = directory.Parent)
            {
                FileInfo candidate = new(Path.Combine(directory.FullName, "SpritzInstaller", "Product.wxs"));
                if (candidate.Exists)
                {
                    return candidate;
                }
            }

            Assert.Fail($"could not find SpritzInstaller/Product.wxs in any ancestor of {startedAt}");
            return null; // unreachable; keeps the compiler happy
        }
    }
}
