using NUnit.Framework;
using RNASeqAnalysisWrappers;
using System.IO;

namespace Test
{
    [TestFixture]
    public class WrapperTests
    {

        #region Installs

        [Test]
        public void test_install_dependencies()
        {
            WrapperUtility.install(TestContext.CurrentContext.TestDirectory);
        }

        [Test]
        public void test_install_star()
        {
            STARWrapper.install(TestContext.CurrentContext.TestDirectory, true, false);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));
        }

        [Test]
        public void test_install_bedops()
        {
            BEDOPSWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));
        }

        [Test]
        public void test_install_rseqc()
        {
            RSeQCWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));
        }

        [Test]
        public void test_install_gatk()
        {
            GATKWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "GenomeAnalysisTK.jar")));
        }

        [Test]
        public void test_install_scalpel()
        {
            ScalpelWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));
        }

        [Test]
        public void test_install_skewer()
        {
            SkewerWrapper.install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "BBMap")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "skewer-0.2.2")));
        }

        #endregion Installs

        #region Minimal alignment tests

        [Test]
        public void test_convert_gff()
        {
            var p = BEDOPSWrapper.gff2bed(Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gff.gff3"), TestContext.CurrentContext.TestDirectory);
            if (p != null) p.WaitForExit();
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample_gff.gff3") + ".bed")).Length > 0);
        }

        [Test]
        public void test_convert_gtf()
        {
            var p = BEDOPSWrapper.gtf2bed(Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"), TestContext.CurrentContext.TestDirectory);
            if (p != null) p.WaitForExit();
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed")).Length > 0);
        }

        [Test]
        public void test_convert_gtf12()
        {
            BEDOPSWrapper.gtf2bed12(Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"), TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test]
        public void test_genome_generate()
        {
            STARWrapper.generate_genome_index(TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.gtf"));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir")));
        }

        [Test]
        public void test_align()
        {
            STARWrapper.basic_align_reads
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

        [Test]
        public void subset_reads_check()
        {
            string[] new_files = new string[0];

            STARWrapper.subset_fastqs(
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

        [Test]
        public void test_strand_specificity_bam()
        {
            Assert.IsTrue(RSeQCWrapper.check_strand_specificity(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.bed"),
                TestContext.CurrentContext.TestDirectory));
        }

        #endregion Minimal alignment tests

        #region Skewer tests

        [Test]
        public void skewer_single()
        {
            SkewerWrapper.trim(TestContext.CurrentContext.TestDirectory, 19,
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
            SkewerWrapper.trim(TestContext.CurrentContext.TestDirectory, 19,
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

        #region STAR-Fusion test

        [Test]
        public void test_star_fusion()
        {
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "fusion_out"));
            STARWrapper.star_fusion(TestContext.CurrentContext.TestDirectory, 
                true, false, 
                8, 
                Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR791578_hg19_Chimeric.out.junction"), 
                new string[0], 
                Path.Combine(TestContext.CurrentContext.TestDirectory, "fusion_out"));
        }

        #endregion STAR-Fusion test

        #region GATK tests

        [Test]
        public void GATK_workflow()
        {
            GATKWrapper.download_known_sites(TestContext.CurrentContext.TestDirectory, TestContext.CurrentContext.TestDirectory, true, true, false, out string known_sites_filename);

            GATKWrapper.prepare_bam(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                out string new_bam);

            GATKWrapper.realign_indels(TestContext.CurrentContext.TestDirectory,
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

            GATKWrapper.variant_calling(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                realigned_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, known_sites_filename),
                out string new_vcf);
        }

        [Test]
        public void variantcall()
        {
            GATKWrapper.variant_calling(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeRep1.Aligned.out.sorted.grouped.marked.split.mapqfixed.realigned.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "common_all_20170710.ensembl.vcf"),
                out string new_vcf);
            Assert.IsTrue(File.Exists(new_vcf));
        }

        #endregion GATK tests

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

    }
}
