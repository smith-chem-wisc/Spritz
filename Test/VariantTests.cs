using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using NUnit.Framework;
using Bio.VCF;
using Genomics;

namespace Test
{
    [TestFixture]
    public class VariantTests
    {
        [Test]
        public void variants_into_genemodel()
        {
            VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.vcf"));
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Assert.AreEqual(15574, variants.Count);

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.gtf"));
            geneModel.amend_transcripts(variants);
        }
    }
}
