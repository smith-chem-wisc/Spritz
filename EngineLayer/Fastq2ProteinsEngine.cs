using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    public class Fastq2ProteinsEngine
    {
        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ucscKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = (int)10e6)
        {
            SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqPaths, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, geneModelGtfOrGff, ucscKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, string[] fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ucscKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = (int)10e6)
        {
            if (Path.GetExtension(genomeFasta) == ".gz")
            {
                WrapperUtility.RunBashCommand("gunzip", WrapperUtility.ConvertWindowsPath(genomeFasta));
                genomeFasta = Path.ChangeExtension(genomeFasta, null);
            }

            string reorderedFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".karyotypic.fa");
            Genome ensemblGenome = new Genome(genomeFasta);
            if (!ensemblGenome.IsKaryotypic())
            {
                if (!File.Exists(reorderedFasta))
                    Genome.WriteFasta(ensemblGenome.KaryotypicOrder(), reorderedFasta);
                ensemblGenome = new Genome(reorderedFasta);
            }
            else
                reorderedFasta = genomeFasta;

            // Parse comma-separated fastq lists
            proteinVariantDatabases = new List<string>();
            if (fastqs.Length > 1 && fastqs[0].Count(x => x == ',') != fastqs[1].Count(x => x == ','))
                return;

            string[] fastqs1 = fastqs[0].Split(',');
            List<string[]> fastqsSeparated = fastqs.Length == 1 ?
                fastqs1.Select(x => new string[] { x }).ToList() :
                fastqs1.Select(x => new string[] { x, fastqs[1].Split(',')[fastqs1.ToList().IndexOf(x)] }).ToList();

            // Trimming
            List<string[]> trimmedFastqSeparated = new List<string[]>();
            List<string> skewerLogs = new List<string>();
            foreach (string[] fq in fastqsSeparated)
            {
                SkewerWrapper.Trim(bin, threads, 19, fq, out string[] trimmedFastqs, out string skewerLog);
                trimmedFastqSeparated.Add(trimmedFastqs);
                skewerLogs.Add(skewerLog);
            }

            // Alignment
            Directory.CreateDirectory(genomeStarIndexDirectory);
            if (!File.Exists(Path.Combine(genomeStarIndexDirectory, "SA")))
            {
                STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff);
            }
            STARWrapper.LoadGenome(bin, genomeStarIndexDirectory);

            List<string> outPrefixes = new List<string>();
            foreach (string[] fq in trimmedFastqSeparated)
            {
                string[] fqForAlignment = fq;
                string outPrefix = Path.Combine(Path.GetDirectoryName(fqForAlignment[0]), Path.GetFileNameWithoutExtension(fqForAlignment[0]));
                outPrefixes.Add(outPrefix);
                if (!File.Exists(outPrefix + STARWrapper.BamFileSuffix) || overwriteStarAlignment)
                {
                    if (inferStrandSpecificity)
                    {
                        STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                        if (useReadSubset) fqForAlignment = subsetFastqs;
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep);
                        strandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                    }
                    STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, fqForAlignment, outPrefix, strandSpecific, STARGenomeLoadOption.LoadAndKeep);
                }
            }
            STARWrapper.LoadGenome(bin, genomeStarIndexDirectory);

            // Variant calling and protein database writing
            List<string> bamPaths = new List<string>();
            Parallel.ForEach(outPrefixes, outPrefix =>
            {
                GATKWrapper.PrepareBamAndFasta(bin, threads, outPrefix + STARWrapper.BamFileSuffix, reorderedFasta, reference, out string newBam);
                lock (bamPaths) bamPaths.Add(newBam);

            });
            foreach (string newBam in bamPaths)
            {
                GATKWrapper.RealignIndels(bin, threads, reorderedFasta, newBam, out string realignedBam);
                GATKWrapper.VariantCalling(bin, threads, reorderedFasta, realignedBam, Path.Combine(bin, ucscKnownSitesPath), out string vcfPath);
                GATKWrapper.ConvertVCFChromosomesUCSC2Ensembl(bin, vcfPath, reference, out string ensemblVcfPath);
                proteinVariantDatabases.Add(
                    WriteSampleSpecificFasta(ensemblVcfPath, ensemblGenome, geneModelGtfOrGff, Path.Combine(Path.GetDirectoryName(fastqs[0]), Path.GetFileNameWithoutExtension(fastqs[0]))));
            }

        }

        public static string WriteSampleSpecificFasta(string vcfPath, Genome genome, string geneModelGtfOrGff, string outprefix)
        {
            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).ToList();
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
            geneModel.AmendTranscripts(singleNucleotideVariants);
            List<Protein> proteins = new List<Protein>();
            Parallel.ForEach(geneModel.Genes, g =>
            {
                List<Protein> someProteins = g.Transcripts.SelectMany(t => t.Translate(true, true)).ToList();
                lock (proteins) proteins.AddRange(someProteins);
            });
            int transcriptsWithVariants = geneModel.Genes.Sum(g => g.Transcripts.Count(x => x.CodingDomainSequences.Any(y => y.Variants.Count > 0)));
            string proteinVariantDatabase =  outprefix + ".protein.fasta";
            ProteinDbWriter.WriteFastaDatabase(proteins, proteinVariantDatabase, " ");
            return proteinVariantDatabase;
        }
    }
}
