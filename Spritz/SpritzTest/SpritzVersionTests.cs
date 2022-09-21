using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;
using System.Linq;

namespace SpritzTest
{
    public class SpritzVersionTests
    {
        [Test]
        public void TestVersionHigher()
        {
            SpritzVersion version = new();
            version.GetVersionNumbersFromWeb();
            bool inReleaseOrlowerVersion =
                version.NewestKnownVersion == RunnerEngine.CurrentVersion && version.NewestKnownVersionWithMsi == null
                || SpritzVersion.IsVersionLower(version.NewestKnownVersionWithMsi);
            Assert.IsTrue(inReleaseOrlowerVersion);
        }

        [Test]
        public void TestVersionMatchInstaller()
        {
            var path = Path.Combine(Environment.CurrentDirectory, @"../../../../../SpritzInstaller/Product.wxs");
            Assert.IsTrue(File.Exists(path));
            var productFile = File.ReadLines(path).ToList();
            string[] linesplit = productFile.First(x => x.Trim().StartsWith("<Product")).Split(' ');
            var version = linesplit.First(x => x.StartsWith("Version")).Split("=")[1].Trim('"');
            bool versionMatchesInstaller = version == RunnerEngine.CurrentVersion;
            Assert.IsTrue(versionMatchesInstaller);
        }
    }
}