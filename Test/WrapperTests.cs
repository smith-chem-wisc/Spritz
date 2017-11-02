using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.IO;
using NUnit.Framework;
using RNASeqAnalysisWrappers;

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

        #endregion Installs

        #region Minimal alignment tests

        [Test]
        public void test_convert_gff()
        {
            var p = BEDOPSWrapper.gff2bed(Path.Combine(TestContext.CurrentContext.TestDirectory, "sample.gff3"), TestContext.CurrentContext.TestDirectory);
            if (p != null) p.WaitForExit();
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample.gff3") + ".bed")).Length > 0);
        }

        [Test]
        public void test_convert_gtf()
        {
            var p = BEDOPSWrapper.gtf2bed(Path.Combine(TestContext.CurrentContext.TestDirectory, "sample.gtf"), TestContext.CurrentContext.TestDirectory);
            if (p != null) p.WaitForExit();
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample.gtf") + ".bed")).Length > 0);
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd1Rep1.fastq.segment.Aligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.bed"),
                TestContext.CurrentContext.TestDirectory));
        }

        #endregion Minimal alignment tests

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
                Path.Combine(TestContext.CurrentContext.TestDirectory,
                "wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd1Rep1.fastq.segment.Aligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                out string new_bam);

            GATKWrapper.realign_indels(TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                new_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, known_sites_filename),
                out string realigned_bam);

            GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa"),
                realigned_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, known_sites_filename),
                out string recal_table_filepath);
        }

        #endregion GATK tests

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
