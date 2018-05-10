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
            ManageToolsFlow.Install(TestContext.CurrentContext.TestDirectory);

            // bedops
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedops")));

            // bedtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "bedtools2", "bin", "bedtools")));

            // cufflinks
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "cufflinks-2.2.1")));

            // gatk
            Assert.IsTrue(Directory.GetFiles(Path.Combine(TestContext.CurrentContext.TestDirectory, "gatk"), "gatk*local.jar").Length > 0);
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "ChromosomeMappings")));

            // hisat2
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "hisat2-2.1.0", "hisat2")));

            // mfold
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "mfold-3.6", "scripts", "mfold")));

            // rsem
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSEM-1.3.0", "rsem-prepare-reference")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSEM-1.3.0", "rsem-calculate-expression")));

            // rseqc
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "RSeQC-2.6.4")));

            // samtools
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "samtools-1.6", "samtools")));

            // scalpel
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "scalpel-0.5.3")));

            // skewer
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "skewer-0.2.2")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "BBMap")));

            // slncky
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky")));
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "slncky", "annotations")));
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "lastz*").Length > 0);

            // snpeff
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "SnpEff", "snpEff.jar")));

            // sratoolkit
            Assert.IsTrue(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*").Length > 0);
            Assert.IsTrue(Directory.GetFiles(Directory.GetDirectories(Directory.GetDirectories(TestContext.CurrentContext.TestDirectory, "sratoolkit*")[0], "bin")[0], "fastq-dump").Length > 0);

            // star
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR")));

            // star-fusion
            Assert.IsTrue(Directory.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "STAR-Fusion_v1.1.0")));
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

        [Test, Order(1)]
        public void TestDownloadSRA()
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            sratoolkit.Fetch(TestContext.CurrentContext.TestDirectory, "SRR6304532", TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(sratoolkit.FastqPaths.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(sratoolkit.LogPath));
        }

        #endregion SRA download test

        #region BED conversion tests

        [Test, Order(1)]
        public void TestConvertGff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void TestConvertGtf12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test, Order(1)]
        public void TestConvertGffToBed12()
        {
            BEDOPSWrapper.Gtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gff.gff3") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        #region Infer Experiment tests

        [Test, Order(4)]
        public void StrandSpecificityTest()
        {
            BAMProperties bam = new BAMProperties(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
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
                new string[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "read1.fastq") },
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
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "read1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","read2.fastq")
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
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "read1.fastq.gz"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "read2.fastq.gz")
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
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                true,
                "grch37");
            Assert.IsTrue(File.Exists(gatk.EnsemblKnownSitesPath));

            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "setupKnownSitesTest.bash");
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
                    Environment.ProcessorCount,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperAgain-trimmedAligned.sortedByCoord.out.bam"),
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

            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "gatkWorkflowTest.bash"), commands).WaitForExit();
            Assert.IsTrue(File.Exists(gatk.HaplotypeCallerVcfPath) && new FileInfo(gatk.HaplotypeCallerVcfPath).Length > 0);
        }

        [Test, Order(5)]
        public void convertVcf()
        {
            var gatk = new GATKWrapper();
            gatk.ConvertVCFChromosomesUCSC2Ensembl(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1Ucsc.vcf"),
                "grch37");
            Assert.IsTrue(File.Exists(gatk.ConvertedEnsemblVcfPath) && new FileInfo(gatk.ConvertedEnsemblVcfPath).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks and Stringtie tests

        [Test, Order(4)]
        public void CufflinksRun()
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "cufflinksRun.bash");
            WrapperUtility.GenerateAndRunScript(script_name, CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
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
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "stringtieRun.bash");
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

        [Test, Order(5)]
        public void SlnckyRun()
        {
            string cufflinksTranscripts = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.cufflinksOutput", CufflinksWrapper.TranscriptsFilename);
            string slnckyOutPrefix = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "cuffmergeModel848835768.slnckyOut", "annotated"); // strange folder, so that it covers the same test as the lncRNAdiscovery run
            string scriptName = Path.Combine(TestContext.CurrentContext.TestDirectory, "scripts", "SlnckyRun.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, 
                SlnckyWrapper.Annotate(
                    TestContext.CurrentContext.TestDirectory,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                    Environment.ProcessorCount,
                    cufflinksTranscripts,
                    "GRCh37",
                    slnckyOutPrefix
                )).WaitForExit();
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(slnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        }

        #endregion Slncky tests

        #region Scalpel tests

        [Test, Order(4)]
        public void ScalpelCall()
        {
            var scalpel = new ScalpelWrapper();
            WrapperUtility.GenerateAndRunScript(Path.Combine(
                TestContext.CurrentContext.TestDirectory, "scripts", "scalpel.bash"),
                scalpel.CallIndels(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.bed12"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "scalpel_test_out")))
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
                "grch37");
            Assert.IsTrue(File.Exists(snpeff.DatabaseListPath));
            string[] databases = Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "snpEff", "data"));
            Assert.IsTrue(databases.Any(x => Path.GetFileName(x).StartsWith("grch37", true, null)));
        }

        [Test, Order(4)]
        public void BasicSnpEffAnnotation()
        {
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "scripts"), 
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    "grch37",
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper-trimmedAligned.sortedByCoord.outProcessed.out.fixedQuals.split.vcf")))
                .WaitForExit();
            Assert.IsTrue(File.Exists(snpeff.HtmlReportPath) && new FileInfo(snpeff.HtmlReportPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedVcfPath) && new FileInfo(snpeff.AnnotatedVcfPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedGenesSummaryPath) && new FileInfo(snpeff.AnnotatedGenesSummaryPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.VariantProteinFastaPath) && new FileInfo(snpeff.VariantProteinFastaPath).Length > 0);
        }

        #endregion SnpEff tests
    }
}