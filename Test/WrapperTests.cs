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
            STARWrapper.install(TestContext.CurrentContext.TestDirectory);
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
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample.gff3") + ".bed")).Length > 0);
        }

        [Test]
        public void test_genome_generate()
        {
            STARWrapper.generate_genome_index(TestContext.CurrentContext.TestDirectory,
                8, 
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1.gtf"));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1")));
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
        public void test_strand_specificity_check()
        {
            //BEDOPSWrapper.gtf2bed("chr1.gtf", TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(STARWrapper.check_strand_specificity(
                TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1"),
                new string[]
                {
                    @"F:\MCF7EncodeData\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd1Rep1.fastq.gz",
                    @"F:\MCF7EncodeData\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd2Rep1.fastq.gz"
                    //Path.Combine(TestContext.CurrentContext.TestDirectory, "r1.fastq"),
                    //Path.Combine(TestContext.CurrentContext.TestDirectory, "r2.fastq")
                },
                "chr1.bed",
                Path.Combine(TestContext.CurrentContext.TestDirectory, TestContext.CurrentContext.TestDirectory)));
        }
    }
}
