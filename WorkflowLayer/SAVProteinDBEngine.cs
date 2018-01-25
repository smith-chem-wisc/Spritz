using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    public class SAVProteinDBEngine
    {

        #region Runner Methods

        public static void GenerateFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            GenerateFromFastqs(bin, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void GenerateFromFastqs(string bin, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            PrepareInputFileFlow.PrepareEnsemblGenomeFasta(genomeFasta, out Genome ensemblGenome, out string reorderedFasta);
            STAR2PassAlignFlow.RunFromFastqs(bin, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, reorderedFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles, out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles, useReadSubset, readSubset);

            proteinVariantDatabases = new List<string>();

            // Variant Calling
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(bin, proteinFasta, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);
            foreach (string dedupedBam in dedupedBamFiles)
            {
                GATKWrapper.VariantCalling(bin, threads, reorderedFasta, dedupedBam, Path.Combine(bin, ensemblKnownSitesPath), out string vcfPath);
                SnpEffWrapper.PrimaryVariantAnnotation(bin, reference, vcfPath, out string htmlReport, out string annotatedVcfPath);
                proteinVariantDatabases.Add(
                    WriteSampleSpecificFasta(annotatedVcfPath, ensemblGenome, badProteinAccessions, selenocysteineContainingAccessions, geneModelGtfOrGff, 7, Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam))));
            }
        }

        #endregion Runner Methods

        #region Sample Specific Database Writing

        public static string WriteSampleSpecificFasta(string vcfPath, Genome genome, HashSet<string> badProteinAccessions, Dictionary<string, string> selenocysteineContaininAccessions, string geneModelGtfOrGff, int minPeptideLength, string outprefix)
        {
            if (!File.Exists(vcfPath))
            {
                Console.WriteLine("Error: VCF not found: " + vcfPath);
                return "Error: VCF not found: " + vcfPath;
            }
            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantSuperContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).Select(v => new VariantSuperContext(v)).ToList();
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
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
                writer.WriteLine(String.Join(",", Transcript.combinatoricFailures) + "\ttranscripts had too many heterozygous variants for combinatorics");
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

    }
}
