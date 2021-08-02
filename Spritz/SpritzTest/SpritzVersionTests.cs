using NUnit.Framework;
using SpritzBackend;
using System;
using System.IO;
using System.Linq;
using YamlDotNet.RepresentationModel;

namespace SpritzTest
{
    public class SpritzVersionTests
    {
        [Test]
        public void TestConfigVersion()
        {
            using var reader = new StreamReader(Path.Combine(Environment.CurrentDirectory, "workflow", "config", "config.yaml"));
            YamlStream yaml = new();
            yaml.Load(reader);
            foreach (var entry in yaml.Documents[0].RootNode as YamlMappingNode)
            {
                if (entry.Key.ToString() == "spritz_version")
                {
                    Assert.AreEqual(RunnerEngine.CurrentVersion, entry.Value.ToString());
                    break;
                }
            }
        }

        [Test]
        public void TestVersionHigher()
        {
            SpritzVersion version = new();
            version.GetVersionNumbersFromWeb();
            Assert.IsTrue(SpritzVersion.IsVersionLower( version.NewestKnownVersionWithMsi));
        }

        [Test]
        public void TestCommonSmkVersion()
        {
            var lines = File.ReadLines(Path.Combine(Environment.CurrentDirectory, "workflow", "rules", "common.smk"));
            var versionLine = lines.FirstOrDefault(line => line.Contains("MIN_SPRITZ_VERSION"));
            Assert.AreEqual(RunnerEngine.CurrentVersion, versionLine.Split('=')[1].Split('#')[0].Trim(new[] { ' ', '"' }));
        }
    }
}