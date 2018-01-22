using Bio;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using WorkflowLayer;

namespace Test
{
    [TestFixture]
    public class WrapperIntegrationTests
    {

        #region Installs

        [Test, Order(0)]
        public void TestInstall()
        {
            InstallFlow.Install(TestContext.CurrentContext.TestDirectory);

            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "snpEff", "snpEff.jar")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR-Fusion_v1.1.0")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));

            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));

            // gatk
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "gatk", "build", "libs", "gatk.jar")));
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

            // mfold
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "mfold-3.6", "scripts", "mfold")));
        }

        private string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa");
        [Test, Order(1)]
        public void DownloadReferences()
        {
            EnsemblDownloadsWrapper.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                out genomeFastaPath,
                out string gtf,
                out string gff,
                out string proteinFasta
            );

            // - a basic set of chromosomes, fairly small ones
            string a = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            // - chromosomes and contigs that test ordering: 9 comes before 22 in karyotipic order, but not lexographic
            string b = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa");

            if (!File.Exists(a) || !File.Exists(b))
            {
                List<ISequence> chromosomes = new FastAParser().Parse(new FileStream(genomeFastaPath, FileMode.Open)).ToList();
                FastAFormatter formatter = new FastAFormatter();

                if (!File.Exists(a))
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("20") || x.ID.StartsWith("21") || x.ID.StartsWith("22")), a);

                if (!File.Exists(b))
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("9") || x.ID.StartsWith("22") || x.ID.StartsWith("GL000210") || x.ID.StartsWith("HG1287_PATCH")), b);
            }

            // Additional setup for small integration tests
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setup.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                "if [ -f Homo_sapiens.GRCh37.75.gtf.gz ]; then gunzip Homo_sapiens.GRCh37.75.gtf.gz; fi",
                @"if [ ! -f 202122.gtf ]; then grep '^20\|^21\|^22' Homo_sapiens.GRCh37.75.gtf > 202122.gtf; fi",
                @"if [ ! -f 922HG1287_PATCH.gtf ]; then grep '^9\|^22\|^HG1287_PATCH\|^GL000210.1' Homo_sapiens.GRCh37.75.gtf > 922HG1287_PATCH.gtf; fi",
            }).WaitForExit();

            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.gtf")));
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
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(2)]
        public void TestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(2)]
        public void TestConvertGtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        #region Minimal STAR alignment tests

        [Test, Order(2)]
        public void SubsetReadsCheck()
        {
            string[] new_files = new string[0];

            STARWrapper.SubsetFastqs(
                TestContext.CurrentContext.TestDirectory,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read2.fastq")
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
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "genomeGenerate.bash"), 
                STARWrapper.GenerateGenomeIndex(TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.gtf"))).WaitForExit();
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir")));
        }

        [Test, Order(3)]
        public void TestAlign()
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "alignReads.bash"),
            STARWrapper.BasicAlignReadCommands
            (
                TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir"),
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read2.fastq")
                },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "r.")
            )).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "r.Aligned.out.bam")));
        }

        #endregion Minimal STAR alignment tests

        #region Tophat alignment tests

        [Test, Order(1)]
        public void TophatAlign()
        {
            TopHatWrapper.GenerateBowtieIndex(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"),
                out string bowtieIndexPrefix);
            Assert.IsTrue(TopHatWrapper.BowtieIndexExists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa")));

            TopHatWrapper.Align(
                TestContext.CurrentContext.TestDirectory,
                bowtieIndexPrefix,
                8,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "mapper.fastq"),
                },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"),
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.bed"),
                0.8));
        }

        [Test, Order(2)]
        public void InnerDistanceTest()
        {
            Assert.AreEqual(132, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "paired_end.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                out string[] outputFiles));
        }

        #endregion Infer Experiment tests

        #region Skewer tests

        [Test, Order(2)]
        public void SkewerSingle()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory,
                1,
                19,
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq") },
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
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read2.fastq")
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

        [Test, Order(2)]
        public void SkewerPairedGz()
        {
            SkewerWrapper.Trim(TestContext.CurrentContext.TestDirectory,
                19,
                1,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "read1.fastq.gz"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","read2.fastq.gz")
                },
                out string[] readTrimmedPaths,
                out string log);
            Assert.True(Path.GetFileName(readTrimmedPaths[0]) == "read1-trimmed-pair1.fastq");
            Assert.True(Path.GetFileName(readTrimmedPaths[1]) == "read1-trimmed-pair2.fastq");
            Assert.True(Path.GetFileName(log) == "read1-trimmed.log");
            File.Delete(readTrimmedPaths[0]);
            File.Delete(readTrimmedPaths[1]);
            File.Delete(log);
        }

        #endregion Skewer tests

        #region GATK tests

        [Test, Order(2)]
        public void DownloadKnownSites()
        {
            GATKWrapper.DownloadEnsemblKnownVariantSites(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                true,
                "grch37",
                out string ensemblKnownSitesPath);
            Assert.IsTrue(File.Exists(ensemblKnownSitesPath));

            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setupKnownSitesTest.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")) + " ]; " +
                    @"then grep '^#\|^chr20\|^chr21\|^chr22\|^20\|^21\|^22' " + WrapperUtility.ConvertWindowsPath(ensemblKnownSitesPath) +
                    " > 202122.vcf; fi",
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")) + " ]; " +
                    @"then grep '^#\|^chr9\|^chr22\|^chrHG1287_PATCH\|chr21_gl000210_random\|^9\|^22\|^HG1287_PATCH\|^GL000210.1' " + WrapperUtility.ConvertWindowsPath(ensemblKnownSitesPath) +
                    " > 922HG1287_PATCH.vcf; fi",
            }).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")));
        }

        [Test, Order(4)]
        public void GatkWorflow()
        {
            GATKWrapper.PrepareBamAndFasta(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                "grch37",
                out string new_bam);

            // No longer needed with HaplotypeCaller
            //GATKWrapper.RealignIndels(TestContext.CurrentContext.TestDirectory,
            //    8,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
            //    new_bam,
            //    out string realigned_bam,
            //    ""); // not including known sites speeds this up substantially, and I'm not planning to use these indels

            // Takes kind of a long time, and it's not recommended for RNA-Seq yet
            //GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122.fa"),
            //    realigned_bam,
            //    out string recal_table_filepath,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122.vcf"));

            GATKWrapper.VariantCalling(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                new_bam,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"),
                out string newVcf);
            Assert.IsTrue(File.Exists(newVcf));
            Assert.IsTrue(new FileInfo(newVcf).Length > 0);
        }

        [Test, Order(5)]
        public void convertVcf()
        {
            GATKWrapper.ConvertVCFChromosomesUCSC2Ensembl(TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1Ucsc.vcf"),
                "grch37",
                out string ensemblVcf);
            Assert.IsTrue(File.Exists(ensemblVcf));
            Assert.IsTrue(new FileInfo(ensemblVcf).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks tests

        [Test, Order(4)]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.out.bam");
            CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                8,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
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

        [Test, Order(4)]
        public void ScalpelCall()
        {
            ScalpelWrapper.CallIndels(TestContext.CurrentContext.TestDirectory,
                8,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.bed12"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "scalpel_test_out"),
                out string newVcf);
            Assert.IsTrue(File.Exists(newVcf));
        }

        #endregion Scalpel tests

        #region SnpEff tests

        [Test, Order(1)]
        public void DownloadSnpEffDatabase()
        {
            SnpEffWrapper.DownloadSnpEffDatabase(TestContext.CurrentContext.TestDirectory,
                "grch37",
                out string databaseListPath);
            Assert.IsTrue(File.Exists(databaseListPath));
            string[] databases = Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "snpEff", "data"));
            Assert.IsTrue(databases.Any(x => Path.GetFileName(x).StartsWith("grch37", true, null)));
        }

        [Test, Order(4)]
        public void BasicSnpEffAnnotation()
        {
            SnpEffWrapper.PrimaryVariantAnnotation(TestContext.CurrentContext.TestDirectory,
                "grch37",
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper-trimmedAligned.sortedByCoord.outProcessed.out.fixedQuals.split.vcf"),
                out string html,
                out string annVcf
                );
            Assert.IsTrue(File.Exists(html) && new FileInfo(html).Length > 0);
            Assert.IsTrue(File.Exists(annVcf) && new FileInfo(annVcf).Length > 0);
        }

        #endregion SnpEff tests

        #region Bigger STAR tests

        //[Test]
        //public void big_genome_generate()
        //{
        //    STARWrapper.generate_genome_index(TestContext.CurrentContext.TestDirectory,
        //        8,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "chr1"),
        //        new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","chr1.fa") },
        //        Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "chr1.gtf"));
        //    Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "chr1")));
        //}

        //[Test]
        //public void big_align()
        //{
        //    STARWrapper.basic_align_reads
        //    (
        //        TestContext.CurrentContext.TestDirectory,
        //        8,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "chr1"),
        //        new string[]
        //        {
        //            @"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\bin\Debug\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd2Rep1.fastq.segment.fastq",
        //            @"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\bin\Debug\wgEncodeCshlLongRnaSeqMcf7CellPapFastqRd1Rep1.fastq.segment.fastq"
        //        },
        //        Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "wgEncodeRep1.")
        //    );
        //    Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "wgEncodeRep1.Aligned.out.bam")));
        //}

        #endregion Bigger STAR tests

        #region Runner Tests

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromFastqs()
        {
            Fastq2ProteinsEngine.RunFromFastqs(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8, 
                new List<string[]>
                {
                    new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapper.fastq") },
                    new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "mapperAgain.fastq") },
                },
                false,
                false,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"),
                out List<string> proteinDatabases);
            foreach (string database in proteinDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));
            }
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromTwoPairsFastqs()
        {
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq"));
            Fastq2ProteinsEngine.RunFromFastqs(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                8,
                new List<string[]>
                {
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000reads_2.fastq") },
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "2000readsAgain_2.fastq") }
                },
                false,
                false,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"),
                out List<string> proteinDatabases);
            foreach (string database in proteinDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                //Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant"))); // no longer happens with variant filtering criteria
            }
        }

        /// <summary>
        /// Handling tough non-karyotypic ordering of chromosomes and an SRA input
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromSRA()
        {
            Fastq2ProteinsEngine.RunFromSra(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                8,
                "SRR6319804",
                false,
                false,
                true,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.gtf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf"), // there is no equivalent of the patch; just checking that that works
                out List<string> proteinDatabases,
                true,
                1000);
            foreach (string database in proteinDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
               // Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant"))); no variants anymore with the filtering criteria
            }
        }

        #endregion Runner Tests

    }
}
