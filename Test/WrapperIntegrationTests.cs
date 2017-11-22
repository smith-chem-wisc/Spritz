using NUnit.Framework;
using RNASeqAnalysisWrappers;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class WrapperIntegrationTests
    {

        #region Installs

        [Test, Order(1)]
        public void TestInstallDependencies()
        {
            WrapperUtility.Install(TestContext.CurrentContext.TestDirectory);
        }

        [Test, Order(2)]
        public void TestInstallSTAR()
        {
            STARWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));
        }

        [Test, Order(2)]
        public void TestInstallSTARFusion()
        {
            STARFusionWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR-Fusion_v1.1.0")));
        }

        [Test, Order(2)]
        public void TestInstallBEDOPS()
        {
            BEDOPSWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));
        }

        [Test, Order(2)]
        public void TestInstallRSeQC()
        {
            RSeQCWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));
        }

        [Test, Order(2)]
        public void TestInstallGATK()
        {
            GATKWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "GenomeAnalysisTK.jar")));
        }

        [Test, Order(2)]
        public void TestInstallScalpel()
        {
            ScalpelWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));
        }

        [Test, Order(2)]
        public void TestInstallSkewer()
        {
            SkewerWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "BBMap")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "skewer-0.2.2")));
        }

        [Test, Order(2)]
        public void TestInstallSRAToolkit()
        {
            SRAToolkitWrapper.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "sratoolkit")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "sratoolkit", "bin", "fastq-dump")));
        }

        #endregion Installs

        #region SRA download test

        [Test, Order(3)]
        public void TestDownloadSRA()
        {
            SRAToolkitWrapper.Fetch(TestContext.CurrentContext.TestDirectory, "SRR6304532", TestContext.CurrentContext.TestDirectory, out string[] fastqs, out string log);
            Assert.IsTrue(fastqs.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(log));
        }

        #endregion SRA download test

        #region Minimal alignment tests

        [Test, Order(3)]
        public void test_convert_gff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(3)]
        public void test_convert_gtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(3)]
        public void test_convert_gtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test, Order(3)]
        public void subset_reads_check()
        {
            string[] new_files = new string[0];

            STARWrapper.SubsetFastqs(
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "r1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "r2.fastq")
                },
                100,
                TestContext.CurrentContext.TestDirectory,
                out new_files);

            foreach (string file in new_files)
            {
                Assert.IsTrue(File.Exists(file));
                Assert.IsTrue(new FileInfo(file).Length > 0);
            }
        }

        [Test, Order(3)]
        public void test_genome_generate()
        {
            STARWrapper.GenerateGenomeIndex(TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.gtf"));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir")));
        }

        [Test, Order(4)]
        public void test_align()
        {
            STARWrapper.BasicAlignReads
            (
                TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir"),
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "r1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "r2.fastq")
                },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "r.")
            );
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "r.Aligned.out.bam")));
        }

        #endregion Minimal alignment tests

        #region Infer Experiment tests

        [Test, Order(3)]
        public void StrandSpecificityTest()
        {
            Assert.IsTrue(RSeQCWrapper.CheckStrandSpecificity(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.bed")));
        }

        [Test, Order(3)]
        public void InnerDistanceTest()
        {
            Assert.AreEqual(1, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.bed")));
        }

        #endregion Infer Experiment tests

        #region Skewer tests

        [Test]
        public void skewer_single()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory, 
                19,
                1,
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "read1.fastq") },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(File.Exists(readTrimmedPaths[0]));
            Assert.True(File.Exists(log));
            File.Delete(readTrimmedPaths[0]);
            File.Delete(log);
        }

        [Test]
        public void skewer_paired()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory, 
                19,
                1,
                new string[] 
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read2.fastq")
                },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(File.Exists(readTrimmedPaths[0]));
            Assert.True(File.Exists(readTrimmedPaths[1]));
            Assert.True(File.Exists(log));
            File.Delete(readTrimmedPaths[0]);
            File.Delete(readTrimmedPaths[1]);
            File.Delete(log);
        }

        #endregion Skewer tests

        #region GATK tests

        [Test]
        public void GATK_workflow()
        {
            GATKWrapper.DownloadKnownSites(TestContext.CurrentContext.TestDirectory, TestContext.CurrentContext.TestDirectory, true, true, false, Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"), out string known_sites_filename);

            GATKWrapper.PrepareBam(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                out string new_bam);

            GATKWrapper.RealignIndels(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                new_bam,
                out string realigned_bam,
                ""); // not including known sites speeds this up substantially

            // Takes kind of a long time, and it's not recommended for RNA-Seq yet
            //GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
            //    realigned_bam,
            //    out string recal_table_filepath,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, known_sites_filename));

            GATKWrapper.VariantCalling(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                realigned_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, known_sites_filename),
                out string new_vcf);
        }

        [Test]
        public void variantcall()
        {
            GATKWrapper.VariantCalling(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "common_all_20170710.ensembl.vcf"),
                out string new_vcf);
            Assert.IsTrue(File.Exists(new_vcf));
        }

        #endregion GATK tests

        #region Cufflinks tests

        [Test]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmedAligned.out.sorted.grouped.marked.split.mapqfixed.bam");
            CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                8,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22.gtf"),
                false,
                true,
                out string outputDirectory
                );
            Assert.IsTrue(File.Exists(Path.Combine(Path.GetDirectoryName(outputDirectory), "transcripts.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(Path.GetDirectoryName(outputDirectory), "skipped.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(Path.GetDirectoryName(outputDirectory), "isoforms.fpkm_tracking")));
            Assert.IsTrue(File.Exists(Path.Combine(Path.GetDirectoryName(outputDirectory), "genes.fpkm_tracking")));
        }

        #endregion Cufflinks tests

        #region Scalpel tests

        [Test]
        public void scalpel_call()
        {
            ScalpelWrapper.call_indels(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.bed"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel_test_out"),
                out string new_vcf);
            Assert.IsTrue(File.Exists(new_vcf));
        }

        #endregion

        #region Bigger STAR tests

        //[Test]
        //public void big_genome_generate()
        //{
        //    STARWrapper.generate_genome_index(TestContext.CurrentContext.TestDirectory,
        //        8,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1"),
        //        new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa") },
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.gtf"));
        //    Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1")));
        //}

        //[Test]
        //public void big_align()
        //{
        //    STARWrapper.basic_align_reads
        //    (
        //        TestContext.CurrentContext.TestDirectory,
        //        8,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1"),
        //        new string[]
        //        {
        //            @"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\bin\Debug\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd2Rep1.fastq.segment.fastq",
        //            @"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\bin\Debug\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd1Rep1.fastq.segment.fastq"
        //        },
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.")
        //    );
        //    Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.bam")));
        //}

        #endregion Bigger STAR tests

        #region Runner Tests

        [Test]
        public void FullProteinRunFromFastqs()
        {
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "setup.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(TestContext.CurrentContext.TestDirectory),
                "if [ ! -f 22.fa ]; then\n  wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.gtf ]; then\n  wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa ]; then\n  gunzip Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.gtf.gz ]; then\n  gunzip Homo_sapiens.GRCh37.75.gtf.gz \nfi",
                "if [ ! -f 22.fa ]; then\n  mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa 22.fa \nfi",
                "if [ ! -f 22.gtf ]; then\n  grep ^22 Homo_sapiens.GRCh37.75.gtf > 22.gtf \nfi"
            }).WaitForExit();

            Fastq2ProteinsRunner.RunFromFastqs(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8, 
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper.fastq") },
                false,
                true,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22.gtf"),
                out string proteinDatabase);
            Assert.IsTrue(new FileInfo(proteinDatabase).Length > 0);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmed.fastq"));
        }

        [Test]
        public void FullProteinRunFromSRA()
        {
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "setup.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(TestContext.CurrentContext.TestDirectory),
                "if [ ! -f 22.fa ]; then\n  wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.gtf ]; then\n  wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa ]; then\n  gunzip Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz \nfi",
                "if [ ! -f Homo_sapiens.GRCh37.75.gtf.gz ]; then\n  gunzip Homo_sapiens.GRCh37.75.gtf.gz \nfi",
                "if [ ! -f 22.fa ]; then\n  mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa 22.fa \nfi",
                "if [ ! -f 22.gtf ]; then\n  grep ^22 Homo_sapiens.GRCh37.75.gtf > 22.gtf \nfi"
            }).WaitForExit();

            Fastq2ProteinsRunner.RunFromSra(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8,
                "SRR6304532",
                false,
                true,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22.gtf"),
                out string proteinDatabase);
            Assert.IsTrue(new FileInfo(proteinDatabase).Length > 0);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmed.fastq"));
        }

        #endregion Runner Tests
    }
}
