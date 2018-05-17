using Bio;
using Bio.IO.FastA;
using NUnit.Framework;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using Bio.VCF;
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
            ManageToolsFlow.Install(TestContext.CurrentContext.TestDirectory);

            // bedops
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "bedops")));

            // bedtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "bedtools2", "bin", "bedtools")));

            // cufflinks
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "cufflinks-2.2.1")));

            // gatk
            Assert.IsTrue(Directory.GetFiles(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "gatk"), "gatk*local.jar").Length > 0);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "ChromosomeMappings")));

            // hisat2
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "hisat2-2.1.0", "hisat2")));

            // mfold
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "mfold-3.6", "scripts", "mfold")));

            // rsem
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "RSEM-1.3.0", "rsem-prepare-reference")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "RSEM-1.3.0", "rsem-calculate-expression")));

            // rseqc
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "RSeQC-2.6.4")));

            // samtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "samtools-1.6", "samtools")));

            // scalpel
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "scalpel-0.5.3")));

            // skewer
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "skewer-0.2.2")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "BBMap")));

            // slncky
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "slncky")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "slncky", "annotations")));
            Assert.IsTrue(Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools"), "lastz*").Length > 0);

            // snpeff
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "SnpEff", "snpEff.jar")));

            // sratoolkit
            Assert.IsTrue(Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools"), "sratoolkit*").Length > 0);
            Assert.IsTrue(Directory.GetFiles(
                Directory.GetDirectories(
                    Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools"), "sratoolkit*")[0], "bin")[0], "fastq-dump").Length > 0);

            // star
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "STAR-" + STARWrapper.STARVersion)));

            // star-fusion
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "STAR-Fusion-v1.4.0", "STAR-Fusion")));

            // trinity
            Assert.IsTrue(Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools"), "trinity*").Length > 0);
        }

        private string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa");

        [Test, Order(1)]
        public void DownloadReferences()
        {
            EnsemblDownloadsWrapper downloadsWrapper = new EnsemblDownloadsWrapper();
            downloadsWrapper.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37");

            // - a basic set of chromosomes, fairly small ones
            string a = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            // - chromosomes and contigs that test ordering: 9 comes before 22 in karyotipic order, but not lexographic
            string b = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa");

            if (!File.Exists(a) || !File.Exists(b))
            {
                List<ISequence> chromosomes = new FastAParser().Parse(new FileStream(genomeFastaPath, FileMode.Open)).ToList();
                FastAFormatter formatter = new FastAFormatter();

                if (!File.Exists(a))
                {
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("20") || x.ID.StartsWith("21") || x.ID.StartsWith("22")), a);
                }
                if (!File.Exists(b))
                {
                    Genome.WriteFasta(chromosomes.Where(x => x.ID.StartsWith("9") || x.ID.StartsWith("22") || x.ID.StartsWith("GL000210") || x.ID.StartsWith("HG1287_PATCH")), b);
                }
            }

            // Additional setup for small integration tests
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Scripts", "setup.bash");
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

        [Test, Order(1)]
        public void TestDownloadSRA()
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            sratoolkit.Fetch(TestContext.CurrentContext.TestDirectory, TestContext.CurrentContext.TestDirectory, "SRR6304532");
            Assert.IsTrue(sratoolkit.FastqPaths.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(sratoolkit.LogPath));
        }

        #endregion SRA download test

        #region BED conversion tests

        [Test, Order(1)]
        public void TestConvertGff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test, Order(1)]
        public void TestConvertGffToBed12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gff.gff3.converted.gff3") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        #region Infer Experiment tests

        [Test, Order(4)]
        public void StrandSpecificityTest()
        {
            BAMProperties bam = new BAMProperties(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                0.8);
            Assert.AreEqual(Strandedness.None, bam.Strandedness);
            Assert.AreEqual(RnaSeqProtocol.SingleEnd, bam.Protocol);
        }

        [Test, Order(2)]
        public void InnerDistanceTest()
        {
            Assert.AreEqual(132, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs"),
                1,
                19,
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read1.fastq") },
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs"),
                19,
                1,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs","read2.fastq")
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs"),
                19,
                1,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read1.fastq.gz"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read2.fastq.gz")
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
            var gatk = new GATKWrapper();
            gatk.DownloadEnsemblKnownVariantSites(
                TestContext.CurrentContext.TestDirectory,
                true,
                "grch37");
            Assert.IsTrue(File.Exists(gatk.EnsemblKnownSitesPath));

            string scriptPath = WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "setupKnownSitesTest.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")) + " ]; " +
                    @"then grep '^#\|^chr20\|^chr21\|^chr22\|^20\|^21\|^22' " + WrapperUtility.ConvertWindowsPath(gatk.EnsemblKnownSitesPath) +
                    " > 202122.vcf; fi",
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")) + " ]; " +
                    @"then grep '^#\|^chr9\|^chr22\|^chrHG1287_PATCH\|chr21_gl000210_random\|^9\|^22\|^HG1287_PATCH\|^GL000210.1' " + WrapperUtility.ConvertWindowsPath(gatk.EnsemblKnownSitesPath) +
                    " > 922HG1287_PATCH.vcf; fi",
            }).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf")));
        }

        [Test, Order(4)]
        public void GatkWorflow()
        {
            var gatk = new GATKWrapper();
            List<string> commands = new List<string>();
            commands.AddRange(gatk.PrepareBamAndFasta(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    Environment.ProcessorCount,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                    "grch37"));
            commands.AddRange(gatk.SplitNCigarReads(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                gatk.PreparedBamPath));

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

            commands.AddRange(gatk.VariantCalling(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                gatk.SplitTrimBamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf")));

            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "gatkWorkflowTest.bash"), commands).WaitForExit();
            Assert.IsTrue(File.Exists(gatk.HaplotypeCallerVcfPath) && new FileInfo(gatk.HaplotypeCallerVcfPath).Length > 0);
        }

        [Test, Order(5)]
        public void convertVcf()
        {
            var gatk = new GATKWrapper();
            var newvcf = gatk.ConvertVCFChromosomesUCSC2Ensembl(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr1Ucsc.vcf"),
                "grch37");
            Assert.IsTrue(File.Exists(newvcf) && new FileInfo(newvcf).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks and Stringtie tests

        [Test, Order(4)]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "cufflinksRun.bash"));
            WrapperUtility.GenerateAndRunScript(script_name, CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                false,
                false,
                out string outputDirectory
                )).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.TranscriptsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.SkippedTranscriptsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.IsoformAbundanceFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(outputDirectory, CufflinksWrapper.GeneAbundanceFilename)));
        }

        [Test, Order(4)]
        public void StringtieRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "stringtieRun.bash"));
            WrapperUtility.GenerateAndRunScript(script_name, StringTieWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa")),
                Strandedness.None,
                false,
                out string outputGtf
                )).WaitForExit();
            Assert.IsTrue(File.Exists(outputGtf));
            Assert.IsTrue(new FileInfo(outputGtf).Length > 0);
        }

        #endregion Cufflinks tests

        #region Slncky tests

        // covered by lncRNAdiscovery workflow test

        //[Test, Order(5)]
        //public void SlnckyRun()
        //{
        //    string cufflinksTranscripts = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.cufflinksOutput", CufflinksWrapper.TranscriptsFilename);
        //    string slnckyOutPrefix = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "cuffmergeModel848835768.slnckyOut", "annotated"); // strange folder, so that it covers the same test as the lncRNAdiscovery run
        //    string scriptName = Path.Combine(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "SlnckyRun.bash"));
        //    WrapperUtility.GenerateAndRunScript(scriptName, 
        //        SlnckyWrapper.Annotate(
        //            TestContext.CurrentContext.TestDirectory,
        //            Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
        //            Environment.ProcessorCount,
        //            cufflinksTranscripts,
        //            "GRCh37",
        //            slnckyOutPrefix
        //        )).WaitForExit();
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
        //    Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        //}

        #endregion Slncky tests

        #region Scalpel tests

        [Test, Order(4)]
        public void ScalpelCall()
        {
            var scalpel = new ScalpelWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "scalpel.bash"),
                scalpel.CallIndels(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.bed12"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "scalpel_test_out")))
            .WaitForExit();
            Assert.IsTrue(File.Exists(scalpel.IndelVcfPath) && new FileInfo(scalpel.IndelVcfPath).Length > 0);
            Assert.IsTrue(File.Exists(scalpel.FilteredIndelVcfPath) && new FileInfo(scalpel.FilteredIndelVcfPath).Length > 0);
        }

        #endregion Scalpel tests

        #region SnpEff tests

        [Test, Order(1)]
        public void DownloadSnpEffDatabase()
        {
            var snpeff = new SnpEffWrapper();
            snpeff.DownloadSnpEffDatabase(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37");
            Assert.IsTrue(File.Exists(snpeff.DatabaseListPath));
            string[] databases = Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "SnpEff", "data"));
            Assert.IsTrue(databases.Any(x => Path.GetFileName(x).StartsWith("grch37", true, null)));
        }

        [Test, Order(4)]
        public void BasicSnpEffAnnotation()
        {
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"), 
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    "grch37",
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper-trimmedAligned.sortedByCoord.outProcessed.out.fixedQuals.split.vcf")))
                .WaitForExit();
            Assert.IsTrue(File.Exists(snpeff.HtmlReportPath) && new FileInfo(snpeff.HtmlReportPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedVcfPath) && new FileInfo(snpeff.AnnotatedVcfPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedGenesSummaryPath) && new FileInfo(snpeff.AnnotatedGenesSummaryPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.VariantProteinFastaPath) && new FileInfo(snpeff.VariantProteinFastaPath).Length > 0);
        }

        #endregion SnpEff tests

        #region Alignment tests
        [Test, Order(1)]
        public void SubsetReadsCheck()
        {
            string[] new_files = new string[0];

            STARWrapper.SubsetFastqs(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs"),
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read2.fastq")
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
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), "SmallGenomeGenerate.bash"),
                STARWrapper.GenerateGenomeIndex(TestContext.CurrentContext.TestDirectory,
                1,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir"),
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa") },
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_one_transcript.gtf"))).WaitForExit();
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir", "SA")));
        }

        [Test, Order(3)]
        public void TestAlign()
        {
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), "SmallAlignReads.bash"),
                STARWrapper.BasicAlignReadCommands
                (
                    TestContext.CurrentContext.TestDirectory,
                    1,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sampleGenomeDir"),
                    new string[]
                    {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "read2.fastq")
                    },
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "r.")
                )).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "r.Aligned.out.bam")));
        }

        [Test, Order(1)]
        public void TophatAlign()
        {
            TopHatWrapper.GenerateBowtieIndex(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"),
                out string bowtieIndexPrefix);
            Assert.IsTrue(TopHatWrapper.BowtieIndexExists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa")));

            TopHatWrapper.Align(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                bowtieIndexPrefix,
                8,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestFastqs", "mapper.fastq"),
                },
                true,
                out string tophatOutDirectory
                );
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatAcceptedHitsFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatAlignmentSummaryFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatDeletionsBEDFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatInsertionsBEDFilename)));
            Assert.IsTrue(File.Exists(Path.Combine(tophatOutDirectory, TopHatWrapper.TophatJunctionsBEDFilename)));
        } 
        #endregion
        #region Workflow Tests

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromFastqs()
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.Reference = "grch37";
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = new List<string[]>
            {
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq") },
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperAgain.fastq") },
            };
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf");
            flow.Parameters.UniProtXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens_202022.xml.gz");
            flow.GenerateSAVProteinsFromFastqs();

            foreach (string database in flow.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains(FunctionalClass.MISSENSE.ToString())));
            }
            foreach (string database in flow.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains(FunctionalClass.MISSENSE.ToString())));
            }
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(4)]
        public void LncRnaDiscoveryRunFromFastqs()
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);
            LncRNADiscoveryFlow flow = new LncRNADiscoveryFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            flow.Parameters.Reference = "grch37";
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = new List<string[]>
            {
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq") },
            };
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.GeneModelGtfOrGff = geneModelPath;
            flow.LncRNADiscoveryFromFastqs();

            Assert.IsTrue(File.Exists(flow.MergedGtfPath));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromTwoPairsFastqs()
        {
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_2.fastq"));

            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.Reference = "grch37";
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = new List<string[]>
            {
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_2.fastq")
                },
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_2.fastq")
                }
            };
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf");
            flow.Parameters.UniProtXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens_202022.xml.gz");

            flow.GenerateSAVProteinsFromFastqs();
            foreach (string database in flow.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                //Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("ANN="))); // no longer see any variations for this test set with variant filtering criteria
            }
            foreach (string database in flow.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                //Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("ANN="))); // no longer see any variations for this test set with variant filtering criteria
            }
        }

        /// <summary>
        /// Handling tough non-karyotypic ordering of chromosomes and an SRA input
        ///
        /// This also tests well-encoded quality scores, so if it starts to fail, check out whether the exit code of the FixMisencodedQualityBaseReads is expected (2 for failure).
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromSRA()
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR6319804");
            var fastqs = SRAToolkitWrapper.GetFastqsFromSras(flow.Parameters.SpritzDirectory, flow.Parameters.AnalysisDirectory, "SRR6319804");
            Directory.CreateDirectory(flow.Parameters.AnalysisDirectory);
            flow.Parameters.Reference = "grch37";
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = fastqs;
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf"); // there is no equivalent of the patch; just checking that that works
            flow.Parameters.UseReadSubset = true;
            flow.Parameters.ReadSubset = 5000;
            flow.GenerateSAVProteinsFromFastqs();

            foreach (string database in flow.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));// no variants anymore with the filtering criteria
            }
            foreach (string database in flow.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));// no variants anymore with the filtering criteria
            }
        } 
        #endregion

        #region RSEM integration tests

        [Test, Order(2)]
        public void RSEMStarCalculate()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.STAR,
                Strandedness.None,
                new[] { newMapper },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        [Test, Order(2)]
        public void RSEMBowtieCalculate()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.Bowtie2,
                Strandedness.None,
                new[] { newMapper },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        [Test, Order(4)]
        public void RSEMStarCalculateFromPaired()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
                RSEMAlignerOption.STAR,
                Strandedness.None,
                new[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR6319804", "SRR6319804_1.segment-trimmed-pair1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR6319804", "SRR6319804_1.segment-trimmed-pair2.fastq"),
                },
                true);
            quantification.QuantifyTranscripts();

            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.IsoformResultsSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GeneResultsSuffix));
            Assert.IsTrue(Directory.Exists(quantification.RsemOutputPrefix + RSEMWrapper.StatDirectorySuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TimeSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.TranscriptBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamSuffix));
            Assert.IsTrue(File.Exists(quantification.RsemOutputPrefix + RSEMWrapper.GenomeSortedBamIndexSuffix));
        }

        // I'm having trouble getting RSEM to work with comma-separated inputs... I think it's because of STAR, which I have had trouble with in this respect in the past.

        //[Test, Order(2)]
        //public void RSEMStarCalculate2Fastq()
        //{
        //    string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar2.fastq");
        //    if (!File.Exists(newMapper))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs",  "mapper.fastq"), newMapper);
        //    }
        //    string newMapper2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar2Again.fastq");
        //    if (!File.Exists(newMapper2))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData",  "TestFastqs","mapper.fastq"), newMapper2);
        //    }

        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
        //        RSEMAlignerOption.STAR,
        //        Strandedness.None,
        //        new[] { newMapper + "," + newMapper2 },
        //        true,
        //        out string referencePrefix,
        //        out string outputPrefix);

        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.IsoformResultsSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GeneResultsSuffix));
        //    Assert.IsTrue(Directory.Exists(outputPrefix + RSEMWrapper.StatDirectorySuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TimeSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamIndexSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamIndexSuffix));
        //}

        //[Test, Order(4)]
        //public void RSEMStarCalculateTwoPairFastq()
        //{
        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf"),
        //        RSEMAlignerOption.STAR,
        //        Strandedness.None,
        //        new[]
        //        {
        //            String.Join(",", new string[] {
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_1.fastq"),
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","TestFastqs", "2000readsAgain_1.fastq") }),
        //            String.Join(",", new string[] {
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","TestFastqs", "2000reads_2.fastq"),
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData","TestFastqs", "2000readsAgain_2.fastq") })
        //        },
        //        true,
        //        out string referencePrefix,
        //        out string outputPrefix);

        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.IsoformResultsSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GeneResultsSuffix));
        //    Assert.IsTrue(Directory.Exists(outputPrefix + RSEMWrapper.StatDirectorySuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TimeSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.TranscriptBamIndexSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamSuffix));
        //    Assert.IsTrue(File.Exists(outputPrefix + RSEMWrapper.GenomeBamIndexSuffix));
        //}

        #endregion

        #region Variant workflow tests

        [Test, Order(3)]
        public void ThereShouldBeUTRs()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "UtrProblem", "gene.gtf"));
            geneModel.ApplyVariants(new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "UtrProblem", "vcf.vcf")).Select(v => new Variant(null, v, genome.Chromosomes[0])).ToList());
            Assert.IsTrue(geneModel.Genes[0].Transcripts[0].UTRs.Count > 0);
        }

        #endregion

        #region Gene Model integration test

        [Test, Order(3)]
        public void FilterTest()
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"));
            EnsemblDownloadsWrapper.FilterGeneModel(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens.GRCh37.75.gtf"),
                genome,
                out string filtered);
            string[] filteredlines = File.ReadAllLines(filtered);
            foreach (string a in filteredlines)
            {
                if (a.StartsWith("#")) { continue; }
                Assert.IsTrue(new[] { "20", "21", "22" }.Contains(a.Split('\t')[0]));
            }
        }

        #endregion
    }
}