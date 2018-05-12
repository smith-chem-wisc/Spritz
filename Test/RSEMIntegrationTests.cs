using NUnit.Framework;
using System;
using System.IO;
using ToolWrapperLayer;
using WorkflowLayer;

namespace Test
{
    [TestFixture]
    public class RSEMIntegrationTests
    {
        [Test, Order(2)]
        public void RSEMStarCalculate()
        {
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperRsemStar.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq"), newMapper);
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
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq"), newMapper);
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
            string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperRsemBowtie.fastq");
            if (!File.Exists(newMapper))
            {
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapper.fastq"), newMapper);
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
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "SRR6319804_1-trimmed-pair1.segment.fastq"),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "SRR6319804_1-trimmed-pair2.segment.fastq"),
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
        //    string newMapper = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperRsemStar2.fastq");
        //    if (!File.Exists(newMapper))
        //    {
        //        File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs",  "mapper.fastq"), newMapper);
        //    }
        //    string newMapper2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "mapperRsemStar2Again.fastq");
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
        //                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TestFastqs", "2000reads_1.fastq"),
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
    }
}