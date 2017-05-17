using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    class TestProteomicsReadWrite
    {
        [Test]
        public void test_read_write_read_xml()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, nice, false, null, out un);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_xml2.xml"), false, nice, false, null, out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));

            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).All(end => end == null || end > 0 && end <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).All(end => end == null || end > 0 && end <= p.Length)));
        }

        [Test]
        public void test_read_Ensembl_pepAllFasta()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>(), ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.xml"));
            List<Protein> ok2 = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.xml"), false, nice, false, null, out un);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));
            Assert.AreEqual("ENSP00000381386", ok[0].Accession);
            Assert.AreEqual("ENSP00000215773", ok[1].Accession);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24316773:-1 gene:ENSG00000099977 transcript:ENST00000398344 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[0].FullName);
            Assert.AreEqual("pep:known chromosome:GRCh37:22:24313554:24322019:-1 gene:ENSG00000099977 transcript:ENST00000350608 gene_biotype:protein_coding transcript_biotype:protein_coding", ok[1].FullName);

            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).All(end => end == null || end > 0 && end <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).All(end => end == null || end > 0 && end <= p.Length)));
        }

        [Test]
        public void test_read_write_read_fasta()
        {
            var nice = new List<Modification>
            {
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null)
            };

            Dictionary<string, Modification> un;
            List<Protein> ok = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression);
            ProteinDbWriter.WriteFastaDatabase(ok, Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), " ");
            List<Protein> ok2 = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, @"rewrite_test_ensembl.pep.all.fasta"), false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression);

            Assert.AreEqual(ok.Count, ok2.Count);
            Assert.True(Enumerable.Range(0, ok.Count).All(i => ok[i].BaseSequence == ok2[i].BaseSequence));

            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).ToArray().All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).ToArray().All(end => end == null || end > 0 && end <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedBeginPosition).ToArray().All(begin => begin == null || begin > 0 && begin <= p.Length)));
            Assert.True(ok2.All(p => p.ProteolysisProducts.Select(x => x.OneBasedEndPosition).ToArray().All(end => end == null || end > 0 && end <= p.Length)));
        }
    }
}
