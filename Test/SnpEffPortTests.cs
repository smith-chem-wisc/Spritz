using Bio;
using Bio.Extensions;
using Bio.VCF;
using NUnit.Framework;
using Proteogenomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class SnpEffPortTests
    {
        [Test]
        public void TestMissenseMutation()
        {
            // Make a transcript
            Sequence seq = new Sequence(Alphabets.DNA, "AAA".Select(cc => (byte)cc).ToArray(), false);
            seq.ID = "1";
            Chromosome c = new Chromosome(seq, null);
            Gene g = new Gene("", c, "+", 1, 3);
            Transcript t = new Transcript("", "", g, "+", 1, 3, "", null);
            Exon x = new Exon(t, seq, 1, 3, seq.ID, "+", null);
            t.Exons = new List<Exon> { x };
            CDS cds = new CDS(t, seq.ID, "+", 1, 3, null, 0);
            t.CodingDomainSequences = new List<CDS> { cds };

            // Make a missense mutation
            // ugh.vcf has a homozygous variation that should change the codon from AAA to AGA, which code for K and R
            // # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
            // 1   2 .   A   G   64.77 . info   GT:AD:DP:GQ:PL  1/1:2,3:5:69:93,0,69
            List<Variant> variants = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "ugh.vcf")).Select(v => new Variant(null, v, new Chromosome(seq, null))).ToList();

            // Make sure it makes it into the DNA sequence
            t.Variants = new HashSet<Variant>(variants);
            List<Transcript> variantTranscripts = GeneModel.ApplyVariantsCombinitorially(t);
            Assert.AreEqual("AAA", SequenceExtensions.ConvertToString(t.Exons[0].Sequence));
            Assert.AreEqual("K", t.Protein().BaseSequence);
            Assert.AreEqual("AGA", SequenceExtensions.ConvertToString(variantTranscripts[0].Exons[0].Sequence));
            Assert.AreEqual("R", variantTranscripts[0].Protein().BaseSequence);

            // Make sure it gets annotated as a missense mutation
            Assert.IsTrue(variantTranscripts[0].VariantAnnotations.Any(str => str.Contains(FunctionalClass.MISSENSE.ToString())));
        }
    }
}