using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;
using System;

namespace WorkflowLayer
{
    public class Fastq2ProteinsEngine
    {

        #region Runner Methods

        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
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

            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "genomeGenerate.bash"), STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff)).WaitForExit();

            List<string[]> fastqsForAlignment = new List<string[]>();
            List<bool> strandSpecificities = new List<bool>();
            List<string> spliceJunctions = new List<string>();
            List<string> sortedBamFiles = new List<string>();
            List<string> dedupedBamFiles = new List<string>();
            proteinVariantDatabases = new List<string>();

            // Trimming and strand specificity
            foreach (string[] fq in fastqs)
            {
                SkewerWrapper.Trim(bin, threads, 19, fq, out string[] trimmedFastqs, out string skewerLog);
                string[] fqForAlignment = trimmedFastqs;

                // Infer strand specificity
                bool localStrandSpecific = strandSpecific;
                if (inferStrandSpecificity || useReadSubset)
                {
                    STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                    if (useReadSubset)
                    {
                        fqForAlignment = subsetFastqs;
                    }
                    if (inferStrandSpecificity)
                    {
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignSubset.bash"), STARWrapper.BasicAlignReadCommands(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep)).WaitForExit();
                        localStrandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                    }
                }
                strandSpecificities.Add(localStrandSpecific);
                fastqsForAlignment.Add(fqForAlignment);
            }

            // Alignment
            List<string> alignmentCommands = new List<string>();
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SpliceJunctionFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.FirstPassAlignmentCommands(bin, threads, genomeStarIndexDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                }
                spliceJunctions.Add(outPrefix + STARWrapper.SpliceJunctionFileSuffix);
            }
            int uniqueSuffix = 1;
            foreach (string f in fastqsForAlignment.SelectMany(f => f))
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, genomeStarIndexDirectory));
            alignmentCommands.AddRange(STARWrapper.ProcessFirstPassSpliceCommands(spliceJunctions, uniqueSuffix, out string spliceJunctionStartDatabase));
            string secondPassGenomeDirectory = genomeStarIndexDirectory + "SecondPass" + uniqueSuffix.ToString();
            alignmentCommands.AddRange(STARWrapper.GenerateGenomeIndex(bin, threads, secondPassGenomeDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff, spliceJunctionStartDatabase));
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SortedBamFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(bin, threads, secondPassGenomeDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));

                }
                sortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                dedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, secondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignReads.bash"), alignmentCommands).WaitForExit();

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
            List<VariantContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).ToList();
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
            geneModel.AmendTranscripts(singleNucleotideVariants);
            List<Protein> proteins = new List<Protein>();
            List<Transcript> transcripts = geneModel.Genes.SelectMany(g => g.Transcripts).ToList();
            for (int i = 0; i < transcripts.Count; i++)
            {
                Console.WriteLine("Processing transcript " + i.ToString() + "/" + transcripts.Count.ToString() + " " + transcripts[i].ID + " " + transcripts[i].ProteinID);
                proteins.AddRange(transcripts[i].Translate(true, true, badProteinAccessions, selenocysteineContaininAccessions));
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
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf(ProteinAnnotation.SingleAminoAcidVariantLabel) > 0).ToString() + "\tSAV sequences");
                List<int> instancesOfSavs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace(ProteinAnnotation.SingleAminoAcidVariantLabel, "").Length) / ProteinAnnotation.SingleAminoAcidVariantLabel.Length).ToList();
                writer.WriteLine(instancesOfSavs.Max().ToString() + "\tmaximum SAVs in a sequence");
                if (instancesOfSavs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSavs.Where(x => x > 0).Average().ToString() + "\taverage SAVs in sequence with one");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf(ProteinAnnotation.SynonymousVariantLabel) > 0).ToString() + "\tsequences with synonymous codon variation");
                List<int> instancesOfSynonymousSnvs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace(ProteinAnnotation.SynonymousVariantLabel, "").Length) / ProteinAnnotation.SingleAminoAcidVariantLabel.Length).ToList();
                writer.WriteLine(instancesOfSynonymousSnvs.Max().ToString() + "\tmaximum synonymous variations in a sequence");
                if (instancesOfSynonymousSnvs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSynonymousSnvs.Where(x => x > 0).Average().ToString() + "\taverage synonymous variations in sequence with one");
            }
            return proteinMetrics;
        }

        #endregion Sample Specific Database Writing

    }
}
