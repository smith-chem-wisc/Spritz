using NUnit.Framework;
using SpritzBackend;
using System.IO;
using System.Linq;

namespace SpritzTest
{
    public class SpritzVersionTests
    {
        // PublishedReleaseIsNotNewerThanCompiledVersion used to live here. It existed to nag about bumping
        // RunnerEngine.CurrentVersion before a release, and it is removed because the version now comes from
        // the git tag, so there is nothing left to forget. Keeping it would have been worse than useless:
        // after releasing 0.3.14 the development default in Directory.Build.props is still 0.3.13, so every
        // subsequent pull request would build at 0.3.13, find a newer published release, and fail until
        // someone bumped the default - reintroducing exactly the manual version chore this change removes.

        /// <summary>
        /// RunnerEngine.CurrentVersion is read from the assembly rather than hardcoded, so this guards the
        /// reflection actually working. If it silently fell back, the GUI would pull a Docker Hub tag like
        /// "smithlab/spritz:0.0.0" - an image that does not exist - and the failure would surface as a
        /// confusing docker pull error rather than as a version problem. Runs offline.
        /// </summary>
        [Test]
        public void CompiledVersionIsReadFromTheAssembly()
        {
            Assert.That(RunnerEngine.CurrentVersion, Does.Match(@"^\d+\.\d+\.\d+"),
                $"CurrentVersion is '{RunnerEngine.CurrentVersion}', which is not a version number; the "
                + "assembly attribute lookup is probably returning nothing");
            Assert.That(RunnerEngine.CurrentVersion, Is.Not.EqualTo("0.0.0"),
                "CurrentVersion fell back to 0.0.0, so neither AssemblyInformationalVersion nor "
                + "AssemblyVersion was readable; check that Directory.Build.props is being imported");
        }

        /// <summary>
        /// The installer's version and RunnerEngine.CurrentVersion must agree, otherwise a release ships an
        /// installer labelled with a different version than the application reports.
        ///
        /// They now agree by construction: both come from the MSBuild Version property, defaulted in
        /// Directory.Build.props and overridden by /p:Version from the git tag at release. So rather than
        /// comparing two literals - which is what this did while both were hand-maintained - this asserts
        /// that the wiring is still in place, i.e. that nobody has hardcoded a version back into Product.wxs.
        /// Runs offline.
        /// </summary>
        [Test]
        public void InstallerVersionComesFromTheBuildRatherThanALiteral()
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
            Assert.That(installerVersion, Is.EqualTo("$(var.Version)"),
                $"{productWxs.Name} declares Version={installerVersion} instead of the WiX preprocessor "
                + "variable $(var.Version). Hardcoding it here reintroduces a second copy of the version that "
                + "has to be bumped by hand and can drift from the git tag and from RunnerEngine.CurrentVersion");
        }

        /// <summary>
        /// The wixproj must actually define the Version preprocessor variable that Product.wxs consumes,
        /// otherwise the installer build fails on an undefined variable. Checked here because the WiX build
        /// itself only runs on Windows, so a break in this wiring would otherwise surface late. Runs offline.
        /// </summary>
        [Test]
        public void InstallerProjectDefinesTheVersionPreprocessorVariable()
        {
            FileInfo productWxs = FindProductWxs();
            FileInfo wixproj = new(Path.Combine(productWxs.DirectoryName, "SpritzInstaller.wixproj"));
            Assert.That(wixproj.Exists, Is.True, $"could not find {wixproj.FullName}");

            string wixprojText = File.ReadAllText(wixproj.FullName);
            Assert.That(wixprojText, Does.Contain("Version=$(Version)"),
                $"{wixproj.Name} does not pass the MSBuild Version property through DefineConstants, so "
                + "$(var.Version) in Product.wxs would be undefined and the installer build would fail");
            Assert.That(wixprojText, Does.Contain("Directory.Build.props"),
                $"{wixproj.Name} no longer imports Directory.Build.props, so it would not see the version "
                + "default and would need its own copy of it");
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
