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
        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqPaths, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, string[] fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            if (Path.GetExtension(genomeFasta) == ".gz")
            {
                WrapperUtility.RunBashCommand("gunzip", WrapperUtility.ConvertWindowsPath(genomeFasta));
                genomeFasta = Path.ChangeExtension(genomeFasta, null);
            }

            // We need to use the same fasta file throughout and have all the VCF and GTF chromosome reference IDs be the same as these.
            // Right now this is based on ensembl references, so those are the chromosome IDs I will be using throughout
            // TODO: try this with UCSC references to judge whether there's a difference in quality / yield / FDR etc in subsequent proteomics analysis
            // This file needs to be in karyotypic order; this allows us not to have to reorder it for GATK analysis
            string ensemblFastaHeaderDelimeter = " ";
            string reorderedFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".karyotypic.fa");
            Genome ensemblGenome = new Genome(genomeFasta);
            if (!ensemblGenome.IsKaryotypic(ensemblFastaHeaderDelimeter))
            {
                ensemblGenome.Chromosomes = ensemblGenome.KaryotypicOrder(ensemblFastaHeaderDelimeter);
                if (!File.Exists(reorderedFasta))
                {
                    Genome.WriteFasta(ensemblGenome.Chromosomes, reorderedFasta);
                }
            }
            else
            {
                reorderedFasta = genomeFasta;
            }

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
                    bool localStrandSpecific = strandSpecific;
                    //if (inferStrandSpecificity)
                    //{
                    //    STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                    //    if (useReadSubset) fqForAlignment = subsetFastqs;
                    //    string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                    //    STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep);
                    //    localStrandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                    //}
                    STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, fqForAlignment, outPrefix, localStrandSpecific, STARGenomeLoadOption.LoadAndKeep);
                }
            }
            STARWrapper.RemoveGenome(bin, genomeStarIndexDirectory);

            // Variant calling and protein database writing
            // todo: check if there are any outPrefixes that are the same and throw an error
            List<string> bamPaths = new List<string>();
            Parallel.ForEach(outPrefixes, outPrefix =>
            {
                GATKWrapper.PrepareBamAndFasta(bin, threads, outPrefix + STARWrapper.BamFileSuffix, reorderedFasta, reference, out string newBam);
                lock (bamPaths) bamPaths.Add(newBam);

            });
            foreach (string newBam in bamPaths)
            {
                GATKWrapper.VariantCalling(bin, threads, reorderedFasta, newBam, Path.Combine(bin, ensemblKnownSitesPath), out string vcfPath);
                proteinVariantDatabases.Add(
                    WriteSampleSpecificFasta(vcfPath, ensemblGenome, geneModelGtfOrGff, Path.Combine(Path.GetDirectoryName(newBam), Path.GetFileNameWithoutExtension(newBam))));
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
