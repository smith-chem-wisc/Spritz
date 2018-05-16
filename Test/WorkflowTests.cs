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

            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.Reference = "grch37";
            flow.Parameters.Threads = Environment.ProcessorCount;
            flow.Parameters.Fastqs = new List<string[]>
            {
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq") },
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperAgain.fastq") },
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
                new[] { Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq") },
            };
            flow.Parameters.GenomeStarIndexDirectory = starIndexDir;
            flow.Parameters.GenomeFasta = genomeFastaPath;
            flow.Parameters.ProteinFasta = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", EnsemblDownloadsWrapper.GRCh37ProteinFastaFilename);
            flow.Parameters.GeneModelGtfOrGff = geneModelPath;
            flow.Parameters.UseReadSubset = true;
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
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_1.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000reads_1.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_1.fastq"));
            if (!File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_2.fastq")))
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000reads_2.fastq"), Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000readsAgain_2.fastq"));

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
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000reads_1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000reads_2.fastq")
                },
                new string[]
                {
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000readsAgain_1.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs","2000readsAgain_2.fastq")
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

            var fastqs = SRAToolkitWrapper.GetFastqsFromSras(TestContext.CurrentContext.TestDirectory, Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData"), "SRR6319804");
            SampleSpecificProteinDBFlow flow = new SampleSpecificProteinDBFlow();
            flow.Parameters.SpritzDirectory = TestContext.CurrentContext.TestDirectory;
            flow.Parameters.AnalysisDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData");
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
    }
}