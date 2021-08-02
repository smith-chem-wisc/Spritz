using NUnit.Framework;
using SpritzModifications;
using System;

namespace SpritzTest
{
    public class SpritzModificationsTests
    {
        [Test]
        public void Test1()
        {
            var mods = ProteinAnnotation.GetUniProtMods(Environment.CurrentDirectory);
            Assert.Greater(mods.Count, 0);
        }
    }
}