using NUnit.Framework;
using SpritzBackend;

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
    }
}
