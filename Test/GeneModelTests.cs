using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using ToolWrapperLayer;

namespace Test
{
    [TestFixture]
    public class GeneModelTests
    {
        Genome genome;

        [OneTimeSetUp]
        public void setup()
        {
            genome = little_genome();
        }

        [Test]
        public void gtfBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.AreEqual(165, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true, false)).ToList();
        }

        [Test]
        public void gffBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.AreEqual(148, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true, false)).ToList();

            //Forward strand, single coding region
            Assert.AreEqual("ENSP00000334393", proteins[0].Accession);
            Assert.AreEqual(
                "MVTEFIFLGLSDSQELQTFLFMLFFVFYGGIVFGNLLIVITVVSDSHLHSPMYFLLANLSLIDLSLSSVTAPKMITDFFSQRKVISFKGCLVQIFLLHFFGGSEMVILIAMGFDRYIAICKPLHYTTIMCGNACVGIMAVTWGIGFLHSVSQLAFAVHLLFCGPNEVDSFYCDLPRVIKLACTDTYRLDIMVIANSGVLTVCSFVLLIISYTIILMTIQHRPLDKSSKALSTLTAHITVVLLFFGPCVFIYAWPFPIKSLDKFLAVFYSVITPLLNPIIYTLRNKDMKTAIRQLRKWDAHSSVKF",
                proteins[0].BaseSequence);

            //Reverse strand, single coding region
            Assert.AreEqual("ENSP00000473460", proteins[1].Accession);
            Assert.AreEqual(
                "TSLWTPQAKLPTFQQLLHTQLLPPSGLFRPSSCFTRAFPGPTFVSWQPSLARFLPVSQQP" +
                "RQAQVLPHTGLSTSSLCLTVASPRPTPVPGHHLRAQNLLKSDSLVPTAASWWPMKAQNLL" +
                "KLTCPGPAPASCQRLQAQPLPHGGFSRPTSSSWLGLQAQLLPHNSLFWPSSCPANGGQCR" +
                "PKTSSSQTLQAHLLLPGGINRPSFDLRTASAGPALASQGLFPGPALASWQLPQAKFLPAC" +
                "QQPQQAQLLPHSGPFRPNL",
                proteins[1].BaseSequence);
        }

        [Test]
        public void gffAppliedToOther()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            GeneModel additional = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_pacbio.gff3"));

            List<Protein> proteins = additional.Genes.SelectMany(g => g.TranslateUsingAnnotatedStartCodons(geneModel, false, 7)).ToList();

            //Forward strand, single coding region
            Assert.AreEqual("PB2015.1.1", proteins[0].Accession);
            Assert.AreEqual(
                "MVTEFIFLGLSDSQELQTFLFMLFFVFYGGIVFGNLLIVITVVSDSHLHSPMYFLLANLSLIDLSLSSVTAPKMITDFFSQRKVISFKGCLVQIFLLHFFGGSEMVILIAMGFDRYIAICKPLHYTTIMCGNACVGIMAVTWGIGFLHSVSQLAFAVHLLFCGPNEVDSFYCDLPRVIKLACTDTYRLDIMVIANSGVLTVCSFVLLIISYTIILMTIQHRPLDKSSKALSTLTAHITVVLLFFGPCVFIYAWPFPIKSLDKFLAVFYSVITPLLNPIIYTLRNKDMKTAIRQLRKWDAHSSVKF",
                proteins[0].BaseSequence);

            //Reverse strand, single coding region
            Assert.AreEqual("PB2015.2.1", proteins[1].Accession);
            Assert.AreEqual(
                "TSLWTPQAKLPTFQQLLHTQLLPPSGLFRPSSCFTRAFPGPTFVSWQPSLARFLPVSQQP" +
                "RQAQVLPHTGLSTSSLCLTVASPRPTPVPGHHLRAQNLLKSDSLVPTAASWWPMKAQNLL" +
                "KLTCPGPAPASCQRLQAQPLPHGGFSRPTSSSWLGLQAQLLPHNSLFWPSSCPANGGQCR" +
                "PKTSSSQTLQAHLLLPGGINRPSFDLRTASAGPALASQGLFPGPALASWQLPQAKFLPAC" +
                "QQPQQAQLLPHSGPFRPNL",
                proteins[1].BaseSequence);
        }

        public Genome little_genome()
        {
            return new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        }

        [Test]
        public void KaryotypicOrder()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headers.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr1" && seqs[1].FriendlyName == "chr2");
            Assert.IsTrue(seqs[22].FriendlyName == "chrX" && seqs[23].FriendlyName == "chrY" && seqs[24].FriendlyName == "chrM");
        }

        [Test]
        public void KaryotypicOrderShort()
        {
            Genome headers = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "headersShort.fa"));
            var seqs = headers.KaryotypicOrder();
            Assert.IsTrue(seqs[0].FriendlyName == "chr9" && seqs[1].FriendlyName == "chr20");
        }

        [Test]
        public void SameProteins()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "chr5script.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ ! -f Homo_sapiens.GRCh38.pep.all.fa ]; then wget ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.pep.all.fa ]; then gunzip " + EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename + ".gz; fi",
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.pep.all.fa"))
            }).WaitForExit();
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename), out var proteinSeqs, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);

            Genome genome = new Genome(@"D:\GRCh38.81\Homo_sapiens.GRCh38.dna.primary_assembly.fa");
            GeneModel geneModel = new GeneModel(genome, @"D:\GRCh38.81\Homo_sapiens.GRCh38.81.gff3");
            List<Protein> geneBasedProteins = geneModel.Translate(true, true, badProteinAccessions, selenocysteineContainingAccessions);
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(@"D:\GRCh38.81\Homo_sapiens.GRCh38.pep.all.fasta", true, DecoyType.None, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression);
            Dictionary<string, string> accSeq = geneBasedProteins.ToDictionary(p => p.Accession, p => p.BaseSequence);

            foreach (Protein p in proteins)
            {
                // now handled with the badAccessions // && !p.BaseSequence.Contains('*') && !seq.Contains('*') && !p.BaseSequence.Contains('X'))
                if (accSeq.TryGetValue(p.Accession, out string seq)) 
                {
                    // there are instances where alternative start codons are used where Ensembl somehow doesn't include them correctly
                    if (!p.FullName.Contains("GRCh38:MT")) Assert.AreEqual(p.BaseSequence, seq);
                    else Assert.AreEqual(p.BaseSequence.Substring(1), seq.Substring(1)); 
                }
            }
        }

        //[Test]
        //public void TryMakingVariantDatabas()
        //{
        //    string gff = @"E:\GffParserCrash\Homo_sapiens.GRCh38.81.gff2.gff3";
        //    string vcf = @"E:\GffParserCrash\978NAT_ACTGAT_L005_R1_001-trimmed-pair1Aligned.out.UCSC.ordered.sorted.grouped.marked.split.mapqfixed.realigned.vcf";
        //    GATKWrapper.ConvertVCFChromosomesUCSC2Ensembl(TestContext.CurrentContext.TestDirectory, vcf, "grch38", out string convertedVcf);
        //    Genome genome = new Genome(@"E:\GffParserCrash\Homo_sapiens.GRCh38.dna.primary_assembly.fa");
        //    genome.Chromosomes = genome.Chromosomes.Where(x => x.ID.Split(' ')[0] == "9").ToList();
        //    Fastq2ProteinsEngine.WriteSampleSpecificFasta(
        //        convertedVcf,
        //        genome,
        //        gff,
        //        Path.Combine(Path.GetDirectoryName(vcf), Path.GetFileNameWithoutExtension(vcf)));
        //}
    }
}
