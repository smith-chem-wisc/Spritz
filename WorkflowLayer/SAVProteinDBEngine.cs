using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    /// <summary>
    /// Create a database of proteins with single amino acid variants (SAVs) using GATK and SnpEff.
    /// </summary>
    public class SAVProteinDBEngine
    {

        #region Runner Methods

        public static void GenerateSAVProteinsFromSra(string binDirectory, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(binDirectory, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            GenerateSAVProteinsFromFastqs(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void GenerateSAVProteinsFromFastqs(string binDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            PrepareEnsemblGenomeFasta(genomeFasta, out Genome ensemblGenome, out string reorderedFasta);
            STAR2PassAlignFlow.AlignFastqs(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, reorderedFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles, out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles, useReadSubset, readSubset);
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(binDirectory, proteinFasta, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);

            // Variant Calling
            string scriptName = Path.Combine(binDirectory, "scripts", "variantCalling.bash");
            List<string> variantCallingCommands = new List<string>();
            List<string> vcfFilePaths = new List<string>();
            List<string> annotatedVcfFilePaths = new List<string>();
            List<string> snpEffHtmlFilePaths = new List<string>();
            List<string> annotatedGenesSummaryPaths = new List<string>();
            foreach (string dedupedBam in dedupedBamFiles)
            {
                variantCallingCommands.AddRange(GATKWrapper.SplitNCigarReads(binDirectory, genomeFasta, dedupedBam, out string splitTrimBam));
                variantCallingCommands.AddRange(GATKWrapper.VariantCalling(binDirectory, threads, reorderedFasta, splitTrimBam, Path.Combine(binDirectory, ensemblKnownSitesPath), out string vcfPath));
                vcfFilePaths.Add(vcfPath);
                variantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(binDirectory, reference, vcfPath, out string htmlReport, out string annotatedVcfPath, out string annotatedGenesSummaryPath));
                annotatedVcfFilePaths.Add(annotatedVcfPath);
                snpEffHtmlFilePaths.Add(htmlReport);
                annotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, variantCallingCommands).WaitForExit();

            // Generate databases
            GeneModel geneModel = new GeneModel(ensemblGenome, geneModelGtfOrGff);
            proteinVariantDatabases = annotatedVcfFilePaths.Select(annotatedVcf => 
                WriteSampleSpecificFasta(annotatedVcf, ensemblGenome, geneModel, badProteinAccessions, selenocysteineContainingAccessions, 7, Path.Combine(Path.GetDirectoryName(annotatedVcf), Path.GetFileNameWithoutExtension(annotatedVcf)))).ToList();
        }

        #endregion Runner Methods

        #region Sample Specific Database Writing

        public static string WriteSampleSpecificFasta(string vcfPath, Genome genome, GeneModel geneModel, HashSet<string> badProteinAccessions, Dictionary<string, string> selenocysteineContaininAccessions, int minPeptideLength, string outprefix)
        {
            if (!File.Exists(vcfPath))
            {
                Console.WriteLine("Error: VCF not found: " + vcfPath);
                return "Error: VCF not found: " + vcfPath;
            }
            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantSuperContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).Select(v => new VariantSuperContext(v)).ToList();
            geneModel.AmendTranscripts(singleNucleotideVariants);
            List<Protein> proteins = new List<Protein>();
            List<Transcript> transcripts = geneModel.Genes.SelectMany(g => g.Transcripts).ToList();
            for (int i = 0; i < transcripts.Count; i++)
            {
                Console.WriteLine("Processing transcript " + i.ToString() + "/" + transcripts.Count.ToString() + " " + transcripts[i].ID + " " + transcripts[i].ProteinID);
                proteins.AddRange(transcripts[i].TranslateFromSnpEffAnnotatedSNVs(true, true, badProteinAccessions, selenocysteineContaininAccessions));
            }
            int transcriptsWithVariants = geneModel.Genes.Sum(g => g.Transcripts.Count(x => x.CodingDomainSequences.Any(y => y.Variants.Count > 0)));
            string proteinVariantDatabase = outprefix + ".protein.fasta";
            string proteinVariantDatabaseXml = outprefix + ".protein.xml";
            List<Protein> proteinsToWrite = proteins.OrderBy(p => p.Accession).Where(p => p.BaseSequence.Length >= minPeptideLength).ToList();
            ProteinDbWriter.WriteFastaDatabase(proteinsToWrite, proteinVariantDatabase, "|");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsToWrite, proteinVariantDatabaseXml);
            WriteProteinFastaMetrics(outprefix, proteinsToWrite);
            return proteinVariantDatabase;
        }

        public static string WriteProteinFastaMetrics(string outprefix, List<Protein> proteinsToWrite)
        {
            string proteinMetrics = outprefix + ".protein.metrics";
            using (StreamWriter writer = new StreamWriter(proteinMetrics))
            {
                writer.WriteLine(Transcript.combinatoricFailures.Count > 0 ? String.Join(",", Transcript.combinatoricFailures) : "0" + "\ttranscripts had too many heterozygous variants for combinatorics");
                Transcript.combinatoricFailures = new List<string>();
                writer.WriteLine(proteinsToWrite.Count.ToString() + "\tprotein sequences");
                writer.WriteLine(proteinsToWrite.Min(p => p.BaseSequence.Length).ToString() + "\tminimum length");
                writer.WriteLine(proteinsToWrite.Max(p => p.BaseSequence.Length).ToString() + "\tmaxium length");
                writer.WriteLine(proteinsToWrite.Average(p => p.BaseSequence.Length).ToString() + "\taverage length");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf("missense") > 0).ToString() + "\tSAV sequences");
                List<int> instancesOfSavs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace("missense", "").Length) / "missense".Length).ToList();
                writer.WriteLine(instancesOfSavs.Max().ToString() + "\tmaximum SAVs in a sequence");
                if (instancesOfSavs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSavs.Where(x => x > 0).Average().ToString() + "\taverage SAVs in sequence with one");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf("synonymous") > 0).ToString() + "\tsequences with synonymous codon variation");
                List<int> instancesOfSynonymousSnvs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace("synonymous", "").Length) / "synonymous".Length).ToList();
                writer.WriteLine(instancesOfSynonymousSnvs.Max().ToString() + "\tmaximum synonymous variations in a sequence");
                if (instancesOfSynonymousSnvs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSynonymousSnvs.Where(x => x > 0).Average().ToString() + "\taverage synonymous variations in sequence with one");
            }
            return proteinMetrics;
        }

        #endregion Sample Specific Database Writing

        #region Preparing Input Files

        /// <summary>
        /// Prepares an Ensembl genome fasta for alignment and all following analysis. The main issue is that Ensembl orders chromosomes lexigraphically, not karyotypically, like some software like GATK expects.
        /// </summary>
        /// <param name="genomeFasta"></param>
        /// <param name="ensemblGenome"></param>
        /// <param name="reorderedFasta"></param>
        private static void PrepareEnsemblGenomeFasta(string genomeFasta, out Genome ensemblGenome, out string reorderedFasta)
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
            reorderedFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".karyotypic.fa");
            ensemblGenome = new Genome(genomeFasta);
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
        }

        #endregion Preparing Input Files

    }
}
