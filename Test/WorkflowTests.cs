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
    public class WorkflowTests
    {
        /// <summary>
        /// Handling multiple fastq files and chromosomes, single-end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromFastqs()
        {
            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);
            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq") },
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperAgain.fastq") },
                },
                false,
                false,
                false,
                starIndexDir,
                genomeFastaPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                geneModelPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"));
            ssdbf.GenerateSAVProteinsFromFastqs();
            foreach (string database in ssdbf.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains(FunctionalClass.MISSENSE.ToString())));
            }
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
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
            LncRNADiscoveryFlow lncRNAdiscovery = new LncRNADiscoveryFlow();
            lncRNAdiscovery.Parameters = new LncRNADiscoveryParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq") },
                },
                false,
                false,
                false,
                starIndexDir,
                genomeFastaPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                geneModelPath,
                true);
            lncRNAdiscovery.LncRNADiscoveryFromFastqs();

            Assert.IsTrue(File.Exists(lncRNAdiscovery.MergedGtfPath));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.CanonicalToLncsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.ClusterInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.FilteredInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.LncsBedSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.LncsInfoSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrfsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrthologsSuffix));
            Assert.IsTrue(File.Exists(lncRNAdiscovery.SlnckyOutPrefix + SlnckyWrapper.OrthologsTopSuffix));
        }

        /// <summary>
        /// Handling multiple fastq files and chromosomes, single end
        /// </summary>
        [Test, Order(3)]
        public void FullProteinRunFromTwoPairsFastqs()
        {
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_2.fastq"));

            string genomeFastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.fa");
            string geneModelPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.gtf");
            string starIndexDir = STARWrapper.GetGenomeStarIndexDirectoryPath(genomeFastaPath, geneModelPath);

            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                TestContext.CurrentContext.TestDirectory,
                "grch37",
                Environment.ProcessorCount,
                new List<string[]>
                {
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000reads_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000reads_2.fastq") },
                    new string[] {
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000readsAgain_1.fastq"),
                        Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000readsAgain_2.fastq") }
                },
                false,
                false,
                false,
                starIndexDir,
                genomeFastaPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                geneModelPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "202122.vcf"));

            ssdbf.GenerateSAVProteinsFromFastqs();
            foreach (string database in ssdbf.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                //Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("ANN="))); // no longer see any variations for this test set with variant filtering criteria
            }
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
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

            var fastqs = SRAToolkitWrapper.GetFastqsFromSras(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), "SRR6319804");
            SampleSpecificProteinDBFlow ssdbf = new SampleSpecificProteinDBFlow();
            ssdbf.Parameters = new SampleSpecificProteinDBParameters(
                TestContext.CurrentContext.TestDirectory,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"),
                "grch37",
                Environment.ProcessorCount,
                fastqs,
                false,
                false,
                false,
                starIndexDir,
                genomeFastaPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename),
                geneModelPath,
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "922HG1287_PATCH.vcf"), // there is no equivalent of the patch; just checking that that works
                null,
                0.05,
                7,
                true,
                5000);
            ssdbf.GenerateSAVProteinsFromFastqs();

            foreach (string database in ssdbf.VariantAnnotatedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));// no variants anymore with the filtering criteria
            }
            foreach (string database in ssdbf.VariantAppliedProteinFastaDatabases)
            {
                Assert.IsTrue(new FileInfo(database).Length > 0);
                Assert.IsTrue(File.ReadAllLines(database).Any(x => x.Contains("variant")));// no variants anymore with the filtering criteria
            }
        }
    }
}