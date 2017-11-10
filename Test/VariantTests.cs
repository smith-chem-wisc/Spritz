using Bio.VCF;
using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class VariantTests
    {
        //[Test]
        //public void variants_into_genemodel()
        //{
        //    VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.vcf"));
        //    List<VariantContext> variants = vcf.Select(x => x).ToList();
        //    Assert.AreEqual(15574, variants.Count);

        //    Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"));
        //    GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.gtf"));
        //    geneModel.amend_transcripts(variants);
        //    List<Protein> proteins = geneModel.genes.SelectMany(g => g.transcripts.SelectMany(t => t.translate(false, true))).ToList();
        //    List<Protein> proteins_without_variants = geneModel.genes.SelectMany(g => g.transcripts.SelectMany(t => t.translate(false, false))).ToList();
        //    HashSet<string> variant_seqs = new HashSet<string>(proteins.Select(p => p.BaseSequence));
        //    HashSet<string> seqs = new HashSet<string>(proteins_without_variants.Select(p => p.BaseSequence));
        //    Assert.IsTrue(variant_seqs.Count != seqs.Count);
        //}

        [Test]
        public void one_transcript_one_homozygous()
        {
            VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr_1_one_homozygous_missense.vcf"));
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Assert.AreEqual(1, variants.Count);

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_one_transcript.gtf"));
            geneModel.AmendTranscripts(variants);
            List<Protein> proteins = geneModel.Translate(true, true).ToList();
            List<Protein> proteins_wo_variant = geneModel.Translate(true, false).ToList();
            Assert.AreEqual(1, geneModel.Genes.Count);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(1, proteins_wo_variant.Count);
            Assert.AreEqual(2, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
            Assert.IsTrue(proteins[0].FullName != null);
            Assert.IsTrue(proteins[0].FullName.Contains(ProteinAnnotation.SingleAminoAcidVariantLabel));
            Assert.IsTrue(proteins[0].FullName.Contains("1:69640"));
        }

        [Test]
        public void one_transcript_one_heterozygous()
        {
            VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr_1_one_heterozygous_missense.vcf"));
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Assert.AreEqual(1, variants.Count);

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_one_transcript.gtf"));
            geneModel.AmendTranscripts(variants);
            List<Protein> proteins = geneModel.Translate(true, true).ToList();
            List<Protein> proteins_wo_variant = geneModel.Translate(true, false).ToList();
            Assert.AreEqual(1, geneModel.Genes.Count);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(1, proteins_wo_variant.Count);
            Assert.AreEqual(2, new HashSet<string> { proteins[0].BaseSequence, proteins[1].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
            Assert.IsTrue(proteins.All(p => p.FullName != null));
            Assert.IsTrue(proteins.Any(p => p.FullName.Contains(ProteinAnnotation.SingleAminoAcidVariantLabel)));
            Assert.IsTrue(proteins.Any(p => p.FullName.Contains("1:69640")));
        }

        [Test]
        public void one_transcript_one_heterozygous_synonymous()
        {
            VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr_1_one_heterozygous_synonymous.vcf"));
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Assert.AreEqual(1, variants.Count);

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_one_transcript.gtf"));
            geneModel.AmendTranscripts(variants);
            List<Protein> proteins = geneModel.Translate(true, true).ToList();
            List<Protein> proteins_wo_variant = geneModel.Translate(true, false).ToList();
            Assert.AreEqual(1, geneModel.Genes.Count);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(1, proteins_wo_variant.Count);
            Assert.AreEqual(1, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
            Assert.IsTrue(!proteins.Any(p => p.FullName.Contains(ProteinAnnotation.SynonymousVariantLabel)));
            Assert.IsTrue(!proteins.Any(p => p.FullName.Contains("1:69666")));
        }

        [Test]
        public void one_transcript_one_homozygous_synonymous()
        {
            VCFParser vcf = new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr_1_one_homozygous_synonymous.vcf"));
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Assert.AreEqual(1, variants.Count);

            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_one_transcript.gtf"));
            geneModel.AmendTranscripts(variants);
            List<Protein> proteins = geneModel.Translate(true, true).ToList();
            List<Protein> proteins_wo_variant = geneModel.Translate(true, false).ToList();
            Assert.AreEqual(1, geneModel.Genes.Count);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(1, proteins_wo_variant.Count);
            Assert.AreEqual(1, new HashSet<string> { proteins[0].BaseSequence, proteins_wo_variant[0].BaseSequence }.Count);
            Assert.IsTrue(proteins.All(p => p.FullName != null));
            Assert.IsTrue(proteins.Any(p => p.FullName.Contains(ProteinAnnotation.SynonymousVariantLabel)));
            Assert.IsTrue(proteins.Any(p => p.FullName.Contains("1:69666")));
        }
    }
}
