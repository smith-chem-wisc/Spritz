using NUnit.Framework;
using WorkflowLayer;
using ToolWrapperLayer;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class WrapperIntegrationTests
    {

        #region Installs

        [Test, Order(0)]
        public void TestInstall()
        {
            InstallFlow.Run(TestContext.CurrentContext.TestDirectory);

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR-Fusion_v1.1.0")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));

            // gatk
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "GenomeAnalysisTK.jar")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "picard.jar")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "ChromosomeMappings")));

            // skewer
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "skewer-0.2.2")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "BBMap")));

            // sratoolkit
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*").Length > 0);
            Assert.IsTrue(Directory.GetFiles(Directory.GetDirectories(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*")[0], "bin")[0], "fastq-dump").Length > 0);

            // slncky
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky", "annotations")));
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "lastz*").Length > 0);
        }

        [Test, Order(1)]
        public void DownloadReferences()
        {
            // Additional setup for small integration tests
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setup.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(TestContext.CurrentContext.TestDirectory),

                "if [ ! -f 22.fa ]; then wget " +
                    "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.20.fa.gz " +
                    "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.21.fa.gz " +
                    "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz " +
                    "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.chromosome.HG1287_PATCH.fa.gz " +
                    "; fi",
                "if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz ]; then gunzip " +
                    "Homo_sapiens.GRCh37.75.dna_sm.chromosome.20.fa.gz " +
                    "Homo_sapiens.GRCh37.75.dna_sm.chromosome.21.fa.gz " +
                    "Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa.gz " +
                    "Homo_sapiens.GRCh37.75.dna_sm.chromosome.HG1287_PATCH.fa.gz " +
                    "; fi",
                "if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.20.fa ]; then mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.20.fa 20.fa; fi",
                "if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.21.fa ]; then mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.21.fa 21.fa; fi",
                "if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa ]; then mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.22.fa 22.fa; fi",
                "if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.HG1287_PATCH.fa ]; then mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.HG1287_PATCH.fa HG1287_PATCH.fa; fi",
                "cat 20.fa 21.fa 22.fa > 202122.fa",
                "cat 22.fa HG1287_PATCH.fa > 22HG1287_PATCH.fa",

                //"if [ ! -f chr1.fa ]; then wget ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz; fi",
                //"if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz ]; then gunzip Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa.gz; fi",
                //"if [ -f Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa ]; then mv Homo_sapiens.GRCh37.75.dna_sm.chromosome.1.fa chr1.fa; fi",

                "if [ ! -f Homo_sapiens.GRCh37.75.gtf ]; then wget ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz; fi",
                "if [ -f Homo_sapiens.GRCh37.75.gtf.gz ]; then gunzip Homo_sapiens.GRCh37.75.gtf.gz; fi",
                //"if [ ! -f chr1.gtf ]; then grep ^1 Homo_sapiens.GRCh37.75.gtf > chr1.gtf; fi",
                @"if [ ! -f 202122.gtf ]; then grep '^20\|^21\|^22' Homo_sapiens.GRCh37.75.gtf > 202122.gtf; fi",
                @"if [ ! -f 22.gtf ]; then grep '^22' Homo_sapiens.GRCh37.75.gtf > 22.gtf; fi",
                @"if [ ! -f 22HG1287_PATCH.gtf ]; then grep '^22\|^HG1287_PATCH' Homo_sapiens.GRCh37.75.gtf > 22HG1287_PATCH.gtf; fi",
            }).WaitForExit();
        }

        #endregion Installs

        #region SRA download test

        [Test, Order(2)]
        public void TestDownloadSRA()
        {
            SRAToolkitWrapper.Fetch(TestContext.CurrentContext.TestDirectory, "SRR6304532", TestContext.CurrentContext.TestDirectory, out string[] fastqs, out string log);
            Assert.IsTrue(fastqs.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(log));
        }

        #endregion SRA download test

        #region BED conversion tests

        [Test, Order(2)]
        public void TestConvertGff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(2)]
        public void TestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(2)]
        public void TestConvertGtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        #region Minimal STAR alignment tests

        [Test, Order(2)]
        public void SubsetReadsCheck()
        {
            string[] new_files = new string[0];

            STARWrapper.SubsetFastqs(
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read2.fastq")
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

        [Test, Order(2)]
        public void TestGenomeGenerate()
        {
            STARWrapper.GenerateGenomeIndex(TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.gtf"));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir")));
        }

        [Test, Order(3)]
        public void TestAlign()
        {
            STARWrapper.BasicAlignReads
            (
                TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "sampleGenomeDir"),
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "read2.fastq")
                },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "r.")
            );
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "r.Aligned.out.bam")));
        }

        #endregion Minimal STAR alignment tests

        #region Tophat alignment tests

        [Test, Order(1)]
        public void TophatAlign()
        {
            TopHatWrapper.GenerateBowtieIndex(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa"),
                out string bowtieIndexPrefix);
            Assert.IsTrue(TopHatWrapper.BowtieIndexExists(Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_sample.fa")));

            TopHatWrapper.Align(
                TestContext.CurrentContext.TestDirectory,
                bowtieIndexPrefix,
                8,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper.fastq"),
                },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "chr1_one_transcript.gtf"),
                true,
                out string tophatOutDirectory
                );
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatAcceptedHitsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatAlignmentSummaryFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatDeletionsBEDFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatInsertionsBEDFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatJunctionsBEDFilename)));
        }

        #endregion Tophat alignment tests

        #region Infer Experiment tests

        [Test, Order(4)]
        public void StrandSpecificityTest()
        {
            Assert.IsFalse(RSeQCWrapper.CheckStrandSpecificity(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmedAligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.bed"),
                0.8));
        }

        [Test, Order(2)]
        public void InnerDistanceTest()
        {
            Assert.AreEqual(132, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "paired_end.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.gtf"),
                out string[] outputFiles));
        }

        #endregion Infer Experiment tests

        #region Skewer tests

        [Test, Order(2)]
        public void SkewerSingle()
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

        [Test, Order(2)]
        public void SkewerPaired()
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

        [Test, Order(1)]
        public void DownloadKnownSites()
        {
            GATKWrapper.DownloadUCSCKnownVariantSites(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                true,
                "grch37",
                out string ucscKnownSitesPath);
            Assert.IsTrue(File.Exists(ucscKnownSitesPath));
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setupKnownSitesTest.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(TestContext.CurrentContext.TestDirectory),
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.vcf")) + " ]; " + 
                    @"then grep '^#\|^chr20\|^chr21\|^chr22' " + WrapperUtility.ConvertWindowsPath(ucscKnownSitesPath) + 
                    " > 202122.vcf; fi",
            }).WaitForExit();
        }

        [Test, Order(3)]
        public void GatkWorflow()
        {
            GATKWrapper.PrepareBamAndFasta(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmedAligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.fa"),
                "grch37",
                out string new_bam,
                out string ucscGenomeFasta);

            GATKWrapper.RealignIndels(TestContext.CurrentContext.TestDirectory,
                8,
                ucscGenomeFasta,
                new_bam,
                out string realigned_bam,
                ""); // not including known sites speeds this up substantially, and I'm not planning to use these indels

            // Takes kind of a long time, and it's not recommended for RNA-Seq yet
            //GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.fa"),
            //    realigned_bam,
            //    out string recal_table_filepath,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.vcf"));

            GATKWrapper.VariantCalling(TestContext.CurrentContext.TestDirectory,
                8,
                ucscGenomeFasta,
                realigned_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.vcf"),
                out string new_vcf);
            Assert.IsTrue(File.Exists(new_vcf));
            Assert.IsTrue(new FileInfo(new_vcf).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks tests

        [Test, Order(4)]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmedAligned.out.bam");
            CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                8,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.gtf"),
                false,
                true,
                out string outputDirectory
                );
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, "transcripts.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, "skipped.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, "isoforms.fpkm_tracking")));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, "genes.fpkm_tracking")));
        }

        #endregion Cufflinks tests

        #region Scalpel tests

        [Test, Order(3)]
        public void ScalpelCall()
        {
            ScalpelWrapper.call_indels(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.bed"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmedAligned.out.bam"),
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

        /// <summary>
        /// Handling multiple chromosomes
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromFastqs()
        {
            Fastq2ProteinsEngine.RunFromFastqs(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8, 
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper.fastq") },
                false,
                true,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.vcf"),
                out string proteinDatabase);
            Assert.IsTrue(new FileInfo(proteinDatabase).Length > 0);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmed.fastq"));
        }

        /// <summary>
        /// Single chromosome, so faster, but 
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromSRA()
        {
            Fastq2ProteinsEngine.RunFromSra(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8,
                "SRR6319804",
                false,
                true,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22HG1287_PATCH"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22HG1287_PATCH.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "22HG1287_PATCH.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "202122.vcf"), // there is no equivalent of the patch; just checking that that works
                out string proteinDatabase,
                true,
                1000);
            Assert.IsTrue(new FileInfo(proteinDatabase).Length > 0);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "mapper-trimmed.fastq"));
        }

        #endregion Runner Tests
    }
}
