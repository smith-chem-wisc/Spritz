using NUnit.Framework;
using SpritzModifications;
using System.IO;

namespace SpritzTest
{
    /// <summary>
    /// Covers the Philosopher/FragPipe header rewrite (issue #221). The header mzLib writes carries the
    /// original raw VCF line as its description, which Philosopher cannot parse; these fix the shape of what
    /// replaces it.
    /// </summary>
    public class FastaHeaderTests
    {
        // A real Spritz variant description, from the mzLib fixture Test/DatabaseTests/HomozygousHLA.xml.
        private const string RawVcfLine =
            @"6\t32637480\t.\tC\tA\t13539.77\t.\tANN=A|missense_variant|MODERATE|HLA-DQA1|ENSG00000196735"
            + @"|transcript|ENST00000395363.5|protein_coding|1/5|c.22C>A|p.Leu8Met|80/1206|22/768|8/255||"
            + @"\tGT:AD:DP:GQ:PL\t1/1:0,337:337:99:13568,1011,0";

        private static string VariantHeader =>
            ">mz|ENST00000395363_L8M|ENST00000395363" + FastaHeaders.VariantTag + RawVcfLine
            + " OS=Homo sapiens GN=HLA-DQA1";

        [Test]
        public void ToPhilosopherHeader_UsesSpPrefixAndDropsTheRawVcfLine()
        {
            string rewritten = FastaHeaders.ToPhilosopherHeader(VariantHeader);

            Assert.That(rewritten, Is.EqualTo(">sp|ENST00000395363_L8M|ENST00000395363 OS=Homo sapiens GN=HLA-DQA1"));
            Assert.That(rewritten, Does.Not.Contain("|missense_variant|"),
                "the raw VCF line must not survive into the header; .vardesc.tsv carries it instead");
            Assert.That(rewritten, Does.Not.Contain(","), "Philosopher mis-parses descriptions containing commas");
        }

        [Test]
        public void ToPhilosopherHeader_KeepsTheVariantAccessionSoTheLookupTableJoins()
        {
            // The accession is the join key into .vardesc.tsv, so the variant suffix has to stay.
            Assert.That(FastaHeaders.ToPhilosopherHeader(VariantHeader), Does.Contain("ENST00000395363_L8M"));
        }

        [Test]
        public void ToPhilosopherHeader_RewritesDecoysToRevSp()
        {
            string decoy = ">rev_mz|ENST00000395363|ENST00000395363  OS=Homo sapiens GN=HLA-DQA1";

            Assert.That(FastaHeaders.ToPhilosopherHeader(decoy),
                Is.EqualTo(">rev_sp|ENST00000395363|ENST00000395363 OS=Homo sapiens GN=HLA-DQA1"),
                "rev_ is Philosopher's default decoy prefix, so the mz->sp swap must not disturb it");
        }

        [Test]
        public void ToPhilosopherHeader_StripsTranscriptPrefix()
        {
            Assert.That(FastaHeaders.ToPhilosopherHeader(">mz|transcript:ENST1|transcript:ENST1 OS=x GN=y"),
                Is.EqualTo(">sp|ENST1|ENST1 OS=x GN=y"));
        }

        [Test]
        public void ToPhilosopherHeader_IsIdempotent()
        {
            string once = FastaHeaders.ToPhilosopherHeader(VariantHeader);

            Assert.That(FastaHeaders.ToPhilosopherHeader(once), Is.EqualTo(once));
        }

        [Test]
        public void ToPhilosopherHeader_LeavesSequenceLinesAlone()
        {
            Assert.That(FastaHeaders.ToPhilosopherHeader("MKVLAAGIVGLNL"), Is.EqualTo("MKVLAAGIVGLNL"));
            Assert.That(FastaHeaders.ToPhilosopherHeader(""), Is.EqualTo(""));
        }

        [Test]
        public void ToPhilosopherHeader_HandlesAHeaderWithNoOrganismField()
        {
            // Truncated annotations happen; the payload then runs to end of line and there is no OS= to stop at.
            Assert.That(FastaHeaders.ToPhilosopherHeader(">mz|ACC_L8M|NAME" + FastaHeaders.VariantTag + RawVcfLine),
                Is.EqualTo(">sp|ACC_L8M|NAME"));
        }

        [Test]
        public void PhilosopherFastaPath_IsASiblingOfTheOriginal()
        {
            string path = FastaHeaders.PhilosopherFastaPath(Path.Combine("out", "combined.spritz.snpeff.protein.fasta"));

            Assert.That(Path.GetFileName(path), Is.EqualTo("combined.spritz.snpeff.protein.fragpipe.fasta"));
            Assert.That(Path.GetDirectoryName(path), Is.EqualTo("out"));
        }

        [Test]
        public void WritePhilosopherFasta_RewritesHeadersAndCopiesSequencesVerbatim()
        {
            string directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "FastaHeaderTests");
            Directory.CreateDirectory(directory);
            string source = Path.Combine(directory, "source.fasta");
            File.WriteAllLines(source, new[] { VariantHeader, "MKVLAAGIVGLNL" });

            string destination = FastaHeaders.PhilosopherFastaPath(source);
            FastaHeaders.WritePhilosopherFasta(source, destination);

            Assert.That(File.ReadAllLines(destination), Is.EqualTo(new[]
            {
                ">sp|ENST00000395363_L8M|ENST00000395363 OS=Homo sapiens GN=HLA-DQA1",
                "MKVLAAGIVGLNL",
            }));
        }
    }
}
