using Bio.VCF;
using NUnit.Framework;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;
using WorkflowLayer;

namespace Test
{
    [TestFixture]
    public class WrapperIntegrationTests
    {
        #region Installs

        [Test, Order(0)]
        public void InstallTest()
        {
            ManageToolsFlow.Install(TestContext.CurrentContext.TestDirectory);
            Assert.IsTrue(WrapperUtility.CheckToolSetup(TestContext.CurrentContext.TestDirectory));
        }

        [Test, Order(1)]
        [TestCase("grch37")]
        [TestCase("grch38")]
        public void EnsemblDownloadReferences(string reference)
        {
            EnsemblDownloadsWrapper downloadsWrapper = new EnsemblDownloadsWrapper();
            downloadsWrapper.DownloadReferences(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                reference);

            // - a basic set of chromosomes, fairly small ones
            string a = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa");

            // - chromosomes and contigs that test ordering: 9 comes before 22 in karyotipic order, but not lexographic
            string b = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + ".fa");

            // Additional setup for small integration tests
            string scriptPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Scripts", "setup.bash");
            if (reference == "grch37")
            {
                WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
                {
                    "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                    SamtoolsWrapper.GetSequencesFromFasta(downloadsWrapper.GenomeFastaPath, new[] { "20", "21", "22" }, a),
                    WrapperUtility.EnsureClosedFileCommands(a),
                    SamtoolsWrapper.GetSequencesFromFasta(downloadsWrapper.GenomeFastaPath, new[] { "9", "22", "GL000210.1"}, b),
                    WrapperUtility.EnsureClosedFileCommands(b),
                    "if [ -f Homo_sapiens.GRCh37.75.gtf.gz ]; then gunzip Homo_sapiens.GRCh37.75.gtf.gz; fi",
                    "if [ ! -f 202122" + reference + @".gtf ]; then grep '^20\|^21\|^22' Homo_sapiens.GRCh37.75.gtf > 202122" + reference + ".gtf; fi",
                    "if [ ! -f 922HG1287_PATCH" + reference + @".gtf ]; then grep '^9\|^22\|^GL000210.1' Homo_sapiens.GRCh37.75.gtf > 922GL" + reference + ".gtf; fi",
                }).WaitForExit();
                Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + ".fa")));
                Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"))));
            }
            else
            {
                WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
                {
                    "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                    SamtoolsWrapper.GetSequencesFromFasta(downloadsWrapper.GenomeFastaPath, new[] { "20", "21", "22" }, a),
                    WrapperUtility.EnsureClosedFileCommands(a),
                    "if [ -f Homo_sapiens.GRCh38.81.gff3.gz ]; then gunzip Homo_sapiens.GRCh38.81.gff3.gz; fi",
                    "if [ ! -f 202122" + reference + @".gff3 ]; then grep '^20\|^21\|^22' Homo_sapiens.GRCh38.81.gff3 > 202122" + reference + ".gff3; fi",
                }).WaitForExit();
            }

            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa")));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"))));
        }

        #endregion Installs

        #region SRA download test

        [Test, Order(1)]
        public void SRAToolkitTestDownload()
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            sratoolkit.Fetch(TestContext.CurrentContext.TestDirectory, TestContext.CurrentContext.TestDirectory, "SRR6304532");
            Assert.IsTrue(sratoolkit.FastqPaths.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(sratoolkit.LogPath));
        }

        #endregion SRA download test

        #region BED conversion tests

        [Test, Order(1)]
        public void BEDOPSTestConvertGff()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void BEDOPSTestConvertGtf()
        {
            string bedPath = BEDOPSWrapper.GtfOrGff2Bed6(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(bedPath).Length > 0);
            File.Delete(bedPath);
        }

        [Test, Order(1)]
        public void BEDOPSTestConvertGtf12()
        {
            BEDOPSWrapper.GffOrGtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gtf.gtf"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gtf.gtf") + ".bed12")).Length > 0);
        }

        [Test, Order(1)]
        public void BEDOPSTestConvertGffToBed12()
        {
            BEDOPSWrapper.GffOrGtf2Bed12(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "sample_gff.gff3"));
            Assert.IsTrue(new FileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", Path.GetFileNameWithoutExtension("sample_gff.gff3.converted.gff3") + ".bed12")).Length > 0);
        }

        #endregion BED conversion tests

        #region Infer Experiment tests

        [Test, Order(2)]
        [TestCase("grch37")]
        public void InnerDistanceTest(string reference)
        {
            Assert.AreEqual(132, RSeQCWrapper.InferInnerDistance(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "paired_end.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", reference.EndsWith("37") ? "202122" + reference + ".gtf" : "202122" + reference + ".gff3"),
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
        [TestCase("grch37")]
        [TestCase("grch38")]
        public void GatkDownloadKnownSites(string reference)
        {
            var gatk = new GATKWrapper();
            gatk.DownloadEnsemblKnownVariantSites(
                TestContext.CurrentContext.TestDirectory,
                true,
                reference);
            Assert.IsTrue(File.Exists(gatk.EnsemblKnownSitesPath));

            string vcf202122Filename = "202122" + reference + ".vcf";
            string vcf922HG1287Filename = "922GL" + reference + ".vcf";
            string scriptPath = WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "setupKnownSitesTest.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData")),
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", vcf202122Filename)) + " ]; " +
                    @"then grep '^#\|^chr20\|^chr21\|^chr22\|^20\|^21\|^22' " + WrapperUtility.ConvertWindowsPath(gatk.EnsemblKnownSitesPath) +
                    " > " + vcf202122Filename + "; fi",
                @"if [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", vcf922HG1287Filename)) + " ]; " +
                    @"then grep '^#\|^chr9\|^chr22\|chr21_gl000210_random\|^9\|^22\|^GL000210.1' " + WrapperUtility.ConvertWindowsPath(gatk.EnsemblKnownSitesPath) +
                    " > " + vcf922HG1287Filename + "; fi",
            }).WaitForExit();
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", vcf202122Filename)));
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", vcf922HG1287Filename)));
        }

        [Test, Order(4)]
        [TestCase("grch37")]
        public void GatkWorflow(string reference)
        {
            var gatk = new GATKWrapper();
            List<string> commands = new List<string>();
            commands.AddRange(gatk.PrepareBamAndFasta(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    Environment.ProcessorCount,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0-trimmedAligned.sortedByCoord.out.bam"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                    reference));
            commands.AddRange(gatk.SplitNCigarReads(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                gatk.PreparedBamPath));

            // No longer needed with HaplotypeCaller
            //GATKWrapper.RealignIndels(TestContext.CurrentContext.TestDirectory,
            //    8,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
            //    new_bam,
            //    out string realigned_bam,
            //    ""); // not including known sites speeds this up substantially, and I'm not planning to use these indels

            // Takes kind of a long time, and it's not recommended for RNA-Seq yet
            //GATKWrapper.base_recalibration(TestContext.CurrentContext.TestDirectory,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122" + reference + ".fa"),
            //    realigned_bam,
            //    out string recal_table_filepath,
            //    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestData", "202122" + reference + ".vcf"));

            commands.AddRange(gatk.VariantCalling(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                gatk.SplitTrimBamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".vcf")));

            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "gatkWorkflowTest.bash"), commands).WaitForExit();
            Assert.IsTrue(File.Exists(gatk.HaplotypeCallerVcfPath) && new FileInfo(gatk.HaplotypeCallerVcfPath).Length > 0);
        }

        [Test, Order(5)]
        [TestCase("grch37")]
        public void GatkConvertVcf(string reference)
        {
            var gatk = new GATKWrapper();
            var newvcf = gatk.ConvertVCFChromosomesUCSC2Ensembl(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "chr1Ucsc.vcf"),
                reference);
            Assert.IsTrue(File.Exists(newvcf) && new FileInfo(newvcf).Length > 0);
        }

        #endregion GATK tests

        #region Cufflinks and Stringtie tests

        [Test, Order(4)]
        [TestCase("grch37")]
        public void CufflinksRun(string reference)
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "cufflinksRun.bash"));
            WrapperUtility.GenerateAndRunScript(script_name, CufflinksWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", reference.EndsWith("37") ? "202122" + reference + ".gtf" : "202122" + reference + ".gff3"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa")),
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
        [TestCase("grch37")]
        public void StringtieRun(string reference)
        {
            string bamPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0-trimmedAligned.sortedByCoord.out.bam");
            string script_name = Path.Combine(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "stringtieRun.bash"));
            WrapperUtility.GenerateAndRunScript(script_name, StringtieWrapper.AssembleTranscripts(
                TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                bamPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", reference.EndsWith("37") ? "202122" + reference + ".gtf" : "202122" + reference + ".gff3"),
                new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa")),
                Strandedness.None,
                false,
                out string outputGtf
                )).WaitForExit();
            Assert.IsTrue(File.Exists(outputGtf));
            Assert.IsTrue(new FileInfo(outputGtf).Length > 0);
        }

        #endregion Cufflinks and Stringtie tests

        #region Slncky tests

        // covered by lncRNAdiscovery workflow test

        //[Test, Order(5)]
        //public void SlnckyRun()
        //{
        //    string cufflinksTranscripts = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0-trimmedAligned.sortedByCoord.out.cufflinksOutput", CufflinksWrapper.TranscriptsFilename);
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
        [TestCase("grch37")]
        public void ScalpelCall(string reference)
        {
            var scalpel = new ScalpelWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "scalpel.bash"),
                scalpel.CallIndels(TestContext.CurrentContext.TestDirectory,
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".bed12"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0-trimmedAligned.sortedByCoord.out.bam"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "scalpel_test_out")))
            .WaitForExit();
            Assert.IsTrue(File.Exists(scalpel.IndelVcfPath) && new FileInfo(scalpel.IndelVcfPath).Length > 0);
            Assert.IsTrue(File.Exists(scalpel.FilteredIndelVcfPath) && new FileInfo(scalpel.FilteredIndelVcfPath).Length > 0);
        }

        #endregion Scalpel tests

        #region SnpEff tests

        [Test, Order(1)]
        [TestCase("grch37")]
        [TestCase("grch38")]
        public void SnpEffDatabaseDownload(string reference)
        {
            var snpeff = new SnpEffWrapper();
            snpeff.DownloadSnpEffDatabase(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                reference);
            Assert.IsTrue(File.Exists(snpeff.DatabaseListPath));
            string[] databases = Directory.GetDirectories(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", "SnpEff", "data"));
            Assert.IsTrue(databases.Any(x => Path.GetFileName(x).StartsWith(reference, true, null)));
        }

        [Test, Order(4)]
        [TestCase("grch37")]
        public void SnpEffAnnotationBasics(string reference)
        {
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"),
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    reference,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs",
                        "mapper0" + reference + "-trimmedAligned.sortedByCoord.outProcessed.out.fixedQuals.split.concat.sorted.vcf"),
                    false))
                .WaitForExit();
            Assert.IsTrue(File.Exists(snpeff.HtmlReportPath) && new FileInfo(snpeff.HtmlReportPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedVcfPath) && new FileInfo(snpeff.AnnotatedVcfPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.AnnotatedGenesSummaryPath) && new FileInfo(snpeff.AnnotatedGenesSummaryPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.VariantProteinFastaPath) && new FileInfo(snpeff.VariantProteinFastaPath).Length > 0);
            Assert.IsTrue(File.Exists(snpeff.VariantProteinXmlPath) && new FileInfo(snpeff.VariantProteinXmlPath).Length > 0);
        }

        [Test, Order(4)]
        [TestCase("grch37", "coding1")]
        [TestCase("grch38", "coding2")]
        public void SnpEffCodingChangeProteins(string reference, string vcfFilename)
        {
            //BuildAndCopySnpeff();
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".snpEffAnnotated.vcf"));
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", vcfFilename + ".vcf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"),
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    reference,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"),
                    true))
                .WaitForExit();

            var fastaLines = File.ReadAllLines(snpeff.VariantProteinFastaPath);
            var xmlProts = ProteinDbLoader.LoadProteinXML(snpeff.VariantProteinXmlPath, true, DecoyType.None, null, false, null, out var un);
            if (reference.EndsWith("37"))
            {
                Assert.IsTrue(fastaLines.Any(l => l.Contains("Val204Ala")));
                Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence == "V" && v.OneBasedBeginPosition == 204 && v.VariantSequence == "A")));
            }
            else
            {
                Assert.IsTrue(fastaLines.Any(l => l.Contains("Leu513Val")));
                Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence == "L" && v.OneBasedBeginPosition == 513 && v.VariantSequence == "V")));
            }
        }

        [Test, Order(4)]
        [TestCase("grch37")]
        public void SnpEffFrameshiftProteins(string reference)
        {
            //BuildAndCopySnpeff();
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, "frameshift1.vcf"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, "frameshift1.snpEffAnnotated.vcf"));
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", "frameshift1.vcf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, "frameshift1.vcf"));
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"),
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    reference,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, "frameshift1.vcf"),
                    true))
                .WaitForExit();
            var fastaLines = File.ReadAllLines(snpeff.VariantProteinFastaPath);
            var xmlProts = ProteinDbLoader.LoadProteinXML(snpeff.VariantProteinXmlPath, true, DecoyType.None, null, false, null, out var un);
            Assert.IsTrue(fastaLines.Any(l => l.Contains("Val204fs")));
            Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.Description.Contains("Val204fs"))));
            Assert.IsTrue(xmlProts.Count(p => p.Accession.Contains("ENST00000316027")) == 1);

            // Frameshift variations should be annotated regarding the protein sequence
            Assert.IsTrue(xmlProts.FirstOrDefault(p => p.SequenceVariations.Any(v => v.Description.Contains("Val204fs")))
                .SequenceVariations.Any(v => v.OriginalSequence.StartsWith("V") && v.VariantSequence.Length > 1));
        }

        [Test, Order(4)]
        [TestCase("grch37", "inframe1")]
        [TestCase("grch38", "inframe2")]
        public void SnpEffInframeInsertionProteins(string reference, string vcfFilename)
        {
            //BuildAndCopySnpeff();
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".snpEffAnnotated.vcf"));
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", vcfFilename + ".vcf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"),
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    reference,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"),
                    true))
                .WaitForExit();
            var xmlProts = ProteinDbLoader.LoadProteinXML(snpeff.VariantProteinXmlPath, true, DecoyType.None, null, false, null, out var un);
            if (reference.EndsWith("37")) { Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence == "V" && v.VariantSequence.Length == 2))); }
            else { Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence.Length == 1 && v.VariantSequence.Length == 2))); }
        }

        [Test, Order(5)] // after lncRNA full run
        [TestCase("grch37", "coding1")]
        [TestCase("grch38", "coding2")]
        public void SnpEffCodingChangeProteinsWithAltGeneModel(string reference, string vcfFilename)
        {
            //BuildAndCopySnpeff();
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".snpEffAnnotated.vcf"));
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", vcfFilename + ".vcf"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"));
            string r = SnpEffWrapper.GenerateDatabase(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                null,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs",
                    reference.EndsWith("37") ? "MergedStringtieModel-806392539.filtered.withcds.gtf" : "MergedStringtieModel-7878119.filtered.withcds.gtf")
                );
            var snpeff = new SnpEffWrapper();
            WrapperUtility.GenerateAndRunScript(
                WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "snpEffTest.bash"),
                snpeff.PrimaryVariantAnnotation(
                    TestContext.CurrentContext.TestDirectory,
                    TestContext.CurrentContext.TestDirectory,
                    r,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestVcfs", reference, vcfFilename + ".vcf"),
                    true))
                .WaitForExit();

            var fastaLines = File.ReadAllLines(snpeff.VariantProteinFastaPath);
            var xmlProts = ProteinDbLoader.LoadProteinXML(snpeff.VariantProteinXmlPath, true, DecoyType.None, null, false, null, out var un);
            if (reference.EndsWith("37"))
            {
                Assert.IsTrue(fastaLines.Any(l => l.Contains("Val204Ala")));
                Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence == "V" && v.OneBasedBeginPosition == 204 && v.VariantSequence == "A")));
            }
            else
            {
                Assert.IsTrue(fastaLines.Any(l => l.Contains("Leu513Val")));
                Assert.IsTrue(xmlProts.Any(p => p.SequenceVariations.Any(v => v.OriginalSequence == "L" && v.OneBasedBeginPosition == 513 && v.VariantSequence == "V")));
            }
        }

        /// <summary>
        /// Automated build and copy from Maven for development of SnpEff extension
        /// </summary>
        private void BuildAndCopySnpeff()
        {
            if (!Directory.Exists(@"C:\Users\Anthony\Documents\GitHub\SnpEff")) { return; }
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(TestContext.CurrentContext.TestDirectory, "SnpEffBuild.bash"),
                new List<string>
                {
                    "export VERSION=4.3; " +
                    "export SNPEFF=/mnt/c/Users/Anthony/Documents/GitHub/SnpEff; " +
                    "sh /mnt/c/Users/Anthony/Documents/GitHub/SnpEff/scripts_build/make.sh; " +
                    "cp /mnt/c/Users/Anthony/Documents/GitHub/SnpEff/target/SnpEff-$VERSION-jar-with-dependencies.jar /mnt/c/Users/Anthony/Documents/GitHub/SnpEff/snpEff.jar; " +
                    "cp /mnt/c/Users/Anthony/Documents/GitHub/SnpEff/snpEff.jar /mnt/e/source/repos/Spritz/Test/bin/Debug/Tools/SnpEff; cd /mnt/e/source/repos/Spritz/Test/bin/Debug/Tools/SnpEff; " +
                    "java -jar  /mnt/e/source/repos/Spritz/Test/bin/Debug/Tools/SnpEff/snpEff.jar"
                }
            ).WaitForExit();
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
        public void STARTestGenomeGenerate()
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
        public void STARTestAlign()
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
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestFastqs", "mapper0.fastq"),
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

        [Test, Order(1)]
        public void Hisat2Align()
        {
            HISAT2Wrapper.GenerateIndex(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa"),
                out string IndexPrefix);
            string genomeFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "chr1_sample.fa");
            Assert.IsTrue(HISAT2Wrapper.IndexExists(genomeFasta));

            HISAT2Wrapper.Align(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                IndexPrefix,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestFastqs", "mapper0.fastq"),
                },
                out string outputDirectory
                );
            var output = outputDirectory;
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", output)));

            HISAT2Wrapper.Align(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                IndexPrefix,
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestFastqs", "mapper0.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory,"TestFastqs", "mapper0.fastq"),
                },
                out string outputDirectoryPaired
                );
            var outputPaired = outputDirectoryPaired;
            Assert.IsTrue(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "Tools", outputPaired)));
        }

        #endregion Alignment tests

        #region Workflow Tests

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(3)]
        [TestCase("grch37", "mapper_chr22indelRegion.fastq")]
        [TestCase("grch38", "mapper_chr22indelRegion.fastq")]
        [TestCase("grch37", "mapper0.fastq")]
        [TestCase("grch38", "mapper0.fastq")]
        [TestCase("grch37", "mapper1.fastq")]
        [TestCase("grch38", "mapper1.fastq")]
        [TestCase("grch37", "mapper2.fastq")]
        [TestCase("grch38", "mapper2.fastq")]
        [TestCase("grch37", "mapper3.fastq")]
        [TestCase("grch38", "mapper3.fastq")]
        public void FullProteinRunFromFastqs(string reference, string fastqFilename)
        {
            BuildAndCopySnpeff();
            string f = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", fastqFilename);
            string fastqPath = Path.Combine(Path.GetDirectoryName(f), Path.GetFileNameWithoutExtension(f) + reference + ".fastq");
            if (!File.Exists(fastqPath)) { File.Copy(f, fastqPath); }
            List<string[]> fastqs = new List<string[]> { new[] { fastqPath } };

            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"));
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs");
            flow.Parameters.Reference = reference;
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = fastqs;
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".vcf");
            flow.Parameters.UniProtXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens_202022.xml.gz");
            flow.Parameters.DoTranscriptIsoformAnalysis = true;
            flow.GenerateSampleSpecificProteinDatabases();

            foreach (string database in flow.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains(FunctionalClass.MISSENSE.ToString())));
            }
            foreach (string database in flow.VariantAnnotatedProteinXmlDatabases)
            {
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(database, true, DecoyType.None, null, false, null, out var un);
                Assert.IsTrue(proteins.Count > 0);
            }
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single end
        /// </summary>
        [Test, Order(3)]
        [TestCase("grch37")]
        public void FullProteinRunFromTwoPairsFastqs(string reference)
        {
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "2000readsAgain_2.fastq"));

            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"));
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.Reference = reference;
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
            flow.Parameters.ProteinFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".vcf");
            flow.Parameters.UniProtXmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "Homo_sapiens_202022.xml.gz");

            flow.GenerateSampleSpecificProteinDatabases();
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
        [TestCase("grch37")]
        public void FullProteinRunFromSRA(string reference)
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + ".fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"));
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "SRR6319804");
            var fastqs = SRAToolkitWrapper.GetFastqsFromSras(flow.Parameters.SpritzDirectory, flow.Parameters.AnalysisDirectory, "SRR6319804");
            Directory.CreateDirectory(flow.Parameters.AnalysisDirectory);
            flow.Parameters.Reference = reference;
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = fastqs;
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.ReferenceGeneModelGtfOrGff = geneModelPath;
            flow.Parameters.EnsemblKnownSitesPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922GL" + reference + ".vcf"); // there is no equivalent of the patch; just checking that that works
            flow.Parameters.UseReadSubset = true;
            flow.Parameters.ReadSubset = 5000;
            flow.GenerateSampleSpecificProteinDatabases();

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

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(4)]
        [TestCase("grch37")]
        public void FullLncRnaDiscoveryRunFromFastqs(string reference)
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3"));
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);
            LncRNADiscoveryFlow flow = new LncRNADiscoveryFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
            flow.Parameters.Reference = reference;
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = new List<string[]>
            {
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0.fastq") },
            };
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.GeneModelGtfOrGff = geneModelPath;
            flow.LncRNADiscoveryFromFastqs();

            Assert.IsTrue(flow.ReconstructedTranscriptModels.All(f => File.Exists(f)));
            Assert.IsTrue(File.Exists(flow.MergedTranscriptModel));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));

            // too small of a test to get these outputs, apparently
            //Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            //Assert.IsTrue(File.Exists(flow.SlnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
        }

        #endregion Workflow Tests

        #region RSEM integration tests

        [Test, Order(2)]
        [TestCase("grch37")]
        public void RSEMStarCalculate(string reference)
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
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
        [TestCase("grch37")]
        public void RSEMStarCalculateGz(string reference)
        {
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperCompressed.fastq"));
            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
                RSEMAlignerOption.STAR,
                Strandedness.None,
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperCompressed.fastq.gz") },
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
        [TestCase("grch37")]
        public void RSEMBowtieCalculate(string reference)
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
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
        [TestCase("grch37")]
        public void RSEMStarCalculateFromPaired(string reference)
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapper0.fastq"), newMapper);
            }

            TranscriptQuantificationFlow quantification = new TranscriptQuantificationFlow();
            quantification.Parameters = new TranscriptQuantificationParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
                Environment.ProcessorCount,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
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
        //[TestCase("grch37")]
        //public void RSEMStarCalculate2Fastq(string reference)
        //{
        //    string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar2.fastq");
        //    if (!File.Exists(newMapper))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs",  "mapper0.fastq"), newMapper);
        //    }
        //    string newMapper2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastqs", "mapperRsemStar2Again.fastq");
        //    if (!File.Exists(newMapper2))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData",  "TestFastqs","mapper0.fastq"), newMapper2);
        //    }

        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
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
        //[TestCase("grch37")]
        //public void RSEMStarCalculateTwoPairFastq(string reference)
        //{
        //    TranscriptQuantificationFlow.QuantifyTranscripts(
        //        TestContext.CurrentContext.TestDirectory,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"),
        //        Environment.ProcessorCount,
        //        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + (reference.EndsWith("37") ? ".gtf" : ".gff3")),
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

        #endregion RSEM integration tests

        #region Variant workflow tests

        [Test, Order(3)]
        [TestCase("grch37")]
        public void VCFParserThereShouldBeUTRs(string reference)
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"));
            GeneModel geneModel = new GeneModel(genome, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "UtrProblem", "gene.gtf"));
            geneModel.ApplyVariants(new VCFParser(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "UtrProblem", "vcf.vcf")).Select(v => new Variant(null, v, genome.Chromosomes[0])).ToList());
            Assert.IsTrue(geneModel.Genes[0].Transcripts[0].UTRs.Count > 0);
        }

        #endregion Variant workflow tests

        #region Gene Model integration test

        [Test, Order(3)]
        [TestCase("grch37")]
        public void EnsemblDownloadsFilterTest(string reference)
        {
            Genome genome = new Genome(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122" + reference + ".fa"));
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

        #endregion Gene Model integration test
    }
}