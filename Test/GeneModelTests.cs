using Bio;
using Bio.Extensions;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public class GeneModelTests
    {
        private Genome genome;

        [OneTimeSetUp]
        public void setup()
        {
            genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
        }

        [Test]
        public void gtfBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.AreEqual(165, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true)).ToList();
        }

        [Test]
        public void gffBasics()
        {
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.AreEqual(148, geneModel.Genes.SelectMany(g => g.Transcripts).Count());
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Translate(true)).ToList();

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
            geneModel.Merge(additional);
            List<Protein> proteins = additional.Genes.SelectMany(g => g.Translate(true)).ToList();

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
        public void TranslateReverseStrand()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript_reverse.gtf"));
            List<Protein> proteins_wo_variant = geneModel.Translate(true).ToList();
            Assert.AreEqual("FFYFIIWSLTLLPRAGLELLTSSDPPASASQSVGITGVSHHAQ",
                proteins_wo_variant[0].BaseSequence);
        }

        [Test, Order(3)]
        public void FilterTest()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"));
            EnsemblDownloadsWrapper.FilterGeneModel(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh37.75.gtf"), genome, out string filtered);
            string[] filteredlines = File.ReadAllLines(filtered);
            foreach (string a in filteredlines)
            {
                if (a.StartsWith("#")) { continue; }
                Assert.IsTrue(new[] { "20", "21", "22" }.Contains(a.Split('\t')[0]));
            }
        }

        [Test]
        public void TranslateAnotherReverseStrand()
        {
            // See http://useast.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000233306;r=7:38362864-38363518;t=ENST00000426402

            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "chr7script.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.7.fa ]; then wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.7.fa ]; then gunzip Homo_sapiens.GRCh38.dna.chromosome.7.fa.gz; fi",
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.7.fa"))
            }).WaitForExit();
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.7.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr7_one_transcript_reverse.gtf"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MQWALAVLLAFLSPASQKSSNLEGRTKSVIRQTGSSAEITCDLAEGSNGYIHWYLHQEGKAPQRLQYYDSYNSKVVLESGVSPGKYYTYASTRNNLRLILRNLIENDFGVYYCATWDG",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateHardReverseStrand()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "chr14script.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.14.fa ]; then wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.14.fa ]; then gunzip Homo_sapiens.GRCh38.dna.chromosome.14.fa.gz; fi",
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.14.fa"))
            }).WaitForExit();
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.14.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "HardReverseStrand", "reverse.gff3"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            ISequence codingSequence = new FastAParser().Parse(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "HardReverseStrand", "codingSeq.fa")).First();
            Assert.AreEqual(SequenceExtensions.ConvertToString(codingSequence),
                SequenceExtensions.ConvertToString(geneModel.Genes[0].Transcripts[0].RetrieveCodingSequence()));
            Assert.AreEqual("MNLQAQPKAQNKRKRCLFGGQEPAPKEQPPPLQPPQQSIRVKEEQYLGHEGPGGAVSTSQ" +
                "PVELPPPSSLALLNSVVYGPERTSAAMLSQQVASVKWPNSVMAPGRGPERGGGGGVSDSS" +
                "WQQQPGQPPPHSTWNCHSLSLYSATKGSPHPGVGVPTYYNHPEALKREKAGGPQLDRYVR" +
                "PMMPQKVQLEVGRPQAPLNSFHAAKKPPNQSLPLQPFQLAFGHQVNRQVFRQGPPPPNPV" +
                "AAFPPQKQQQQQQPQQQQQQQQAALPQMPLFENFYSMPQQPSQQPQDFGLQPAGPLGQSH" +
                "LAHHSMAPYPFPPNPDMNPELRKALLQDSAPQPALPQVQIPFPRRSRRLSKEGILPPSAL" +
                "DGAGTQPGQEATGNLFLHHWPLQQPPPGSLGQPHPEALGFPLELRESQLLPDGERLAPNG" +
                "REREAPAMGSEEGMRAVSTGDCGQVLRGGVIQSTRRRRRASQEANLLTLAQKAVELASLQ" +
                "NAKDGSGSEEKRKSVLASTTKCGVEFSEPSLATKRAREDSGMVPLIIPVSVPVRTVDPTE" +
                "AAQAGGLDEDGKGPEQNPAEHKPSVIVTRRRSTRIPGTDAQAQAEDMNVKLEGEPSVRKP" +
                "KQRPRPEPLIIPTKAGTFIAPPVYSNITPYQSHLRSPVRLADHPSERSFELPPYTPPPIL" +
                "SPVREGSGLYFNAIISTSTIPAPPPITPKSAHRTLLRTNSAEVTPPVLSVMGEATPVSIE" +
                "PRINVGSRFQAEIPLMRDRALAAADPHKADLVWQPWEDLESSREKQRQVEDLLTAACSSI" +
                "FPGAGTNQELALHCLHESRGDILETLNKLLLKKPLRPHNHPLATYHYTGSDQWKMAERKL" +
                "FNKGIAIYKKDFFLVQKLIQTKTVAQCVEFYYTYKKQVKIGRNGTLTFGDVDTSDEKSAQ" +
                "EEVEVDIKTSQKFPRVPLPRRESPSEERLEPKREVKEPRKEGEEEVPEIQEKEEQEEGRE" +
                "RSRRAAAVKATQTLQANESASDILILRSHESNAPGSAGGQASEKPREGTGKSRRALPFSE" +
                "KKKKTETFSKTQNQENTFPCKKCGR",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateSelenocysteineContaining()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "chr5script.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ ! -f Homo_sapiens.GRCh38.pep.all.fa ]; then wget ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.pep.all.fa ]; then gunzip " + EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename + ".gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.5.fa ]; then wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.5.fa ]; then gunzip Homo_sapiens.GRCh38.dna.chromosome.5.fa.gz; fi",
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.5.fa")),
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.pep.all.fa"))
            }).WaitForExit();
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.5.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr5_selenocysteineContaining.gff3"));
            var d = new EnsemblDownloadsWrapper();
            d.GetImportantProteinAccessions(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename));
            List<Protein> proteins = geneModel.Translate(true, d.BadProteinAccessions, d.SelenocysteineProteinAccessions).ToList();
            Assert.AreEqual("MWRSLGLALALCLLPSGGTESQDQSSLCKQPPAWSIRDQDPMLNSNGSVTVVALLQASUYLCILQASKLEDLRVKLKKEGYSNISYIVVNHQGISSRLKYTHLKNKVSEHIPVYQQEENQTDVWTLLNGSKDDFLIYDRCGRLVYHLGLPFSFLTFPYVEEAIKIAYCEKKCGNCSLTTLKDEDFCKRVSLATVDKTVETPSPHYHHEHHHNHGHQHLGSSELSENQQPGAPNAPTHPAPPGLHHHHKHKGQHRQGHPENRDMPASEDLQDLQKKLCRKRCINQLLCKLPTDSELAPRSUCCHCRHLIFEKTGSAITUQCKENLPSLCSUQGLRAEENITESCQURLPPAAUQISQQLIPTEASASURUKNQAKKUEUPSN",
                proteins[0].BaseSequence);
        }

        [Test]
        public void TranslateMTSeq()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "chr7script.bash"), new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.MT.fa ]; then wget ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz; fi",
                "if [ ! -f Homo_sapiens.GRCh38.dna.chromosome.MT.fa ]; then gunzip Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz; fi",
                WrapperUtility.EnsureClosedFileCommands(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.MT.fa"))
            }).WaitForExit();
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh38.dna.chromosome.MT.fa"));

            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chrM_one_transcript_reverse.gtf"));
            List<Protein> proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MPMANLLLLIVPILIAMAFLMLTERKILGYMQLRKGPNVVGPYGLLQPFADAMKLFTKEPLKPATSTITLYITAPTLALTIALLLWTPLPMPNPLVNLNLGLLFILATSSLAVYSILWSGWASNSNYALIGALRAVAQTISYE" + 
                "VTLAIILLSTLLMSGSFNLSTLITTQEHLWLLLPSWPLAMMWFISTLAETNRTPFDLAEGESELVSGFNIEYAAGPFALFFMAEYTNIIMMNTLTTTIFLGTTYDALSPELYTTYFVTKTLLLTSLFLWIRTAYPRFRYDQLMHLLWKNFLPLTLALLMWYVSMPITISSIPPQT",
                proteins[0].BaseSequence);

            geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chrM_one_transcript_reverse2.gtf"));
            proteins = geneModel.Translate(true).ToList();
            Assert.AreEqual("MNPLAQPVIYSTIFAGTLITALSSHWFFTWVGLEMNMLAFIPVLTKKMNPRSTEAAIKYFLTQATASMILLMAILFNNMLSGQWTMTNTTNQYSSLMIMMAMAMKLGMAPFHFWVPEVTQGTPLTSGLLLLTWQKLAPISIMYQISPS" + 
                "LNVSLLLTLSILSIMAGSWGGLNQTQLRKILAYSSITHMGWMMAVLPYNPNMTILNLTIYIILTTTAFLLLNLNSSTTTLLLSRTWNKLTWLTPLIPSTLLSLGGLPPLTGFLPKWAIIEEFTKNNSLIIPTIMATITLLNLYFYLRLIYSTSITLLPMSNNVKM" + 
                "KWQFEHTKPTPFLPTLIALTTLLLPISPFMLMIL",
                proteins[0].BaseSequence);
        }

        [Test]
        public void SameProteins()
        {
            EnsemblDownloadsWrapper downloads = new EnsemblDownloadsWrapper();
            downloads.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch38");
            downloads.GetImportantProteinAccessions(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh38ProteinFastaFilename));

            Genome genome = new Genome(downloads.GenomeFastaPath);
            GeneModel geneModel = new GeneModel(genome, downloads.Gff3GeneModelPath);
            List<Protein> geneBasedProteins = geneModel.Translate(true, downloads.BadProteinAccessions, downloads.SelenocysteineProteinAccessions);
            List<Protein> pepDotAll = ProteinDbLoader.LoadProteinFasta(downloads.ProteinFastaPath, true, DecoyType.None, false, 
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out List<string> errors);
            Dictionary<string, string> accSeq = geneBasedProteins.ToDictionary(p => p.Accession, p => p.BaseSequence);

            foreach (Protein p in pepDotAll)
            {
                // now handled with the badAccessions // && !p.BaseSequence.Contains('*') && !seq.Contains('*') && !p.BaseSequence.Contains('X'))
                if (accSeq.TryGetValue(p.Accession, out string seq))
                {
                    Assert.AreEqual(p.BaseSequence, seq);
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