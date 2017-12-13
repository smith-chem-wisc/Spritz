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
        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases, bool useReadSubset = false, int readSubset = 300000)
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
            Directory.CreateDirectory(genomeStarIndexDirectory);
            if (!File.Exists(Path.Combine(genomeStarIndexDirectory, "SA")))
            {
                STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff);
            }

            proteinVariantDatabases = new List<string>();
            STARWrapper.LoadGenome(bin, genomeStarIndexDirectory);
            foreach (string[] fq in fastqs)
            {
                // Trimming
                SkewerWrapper.Trim(bin, threads, 19, fq, out string[] trimmedFastqs, out string skewerLog);

                // Alignment
                string[] fqForAlignment = trimmedFastqs;
                string outPrefix = Path.Combine(Path.GetDirectoryName(fqForAlignment[0]), Path.GetFileNameWithoutExtension(fqForAlignment[0]));
                if (!File.Exists(outPrefix + STARWrapper.BamFileSuffix) || overwriteStarAlignment)
                {
                    bool localStrandSpecific = strandSpecific;
                    if (inferStrandSpecificity)
                    {
                        STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                        if (useReadSubset) fqForAlignment = subsetFastqs;
                        string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                        STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep);
                        localStrandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                    }
                    STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, fqForAlignment, outPrefix, localStrandSpecific, STARGenomeLoadOption.LoadAndKeep);
                }

                // Variant Calling
                GATKWrapper.PrepareBamAndFasta(bin, threads, outPrefix + STARWrapper.BamFileSuffix, reorderedFasta, reference, out string newBam);
                GATKWrapper.VariantCalling(bin, threads, reorderedFasta, newBam, Path.Combine(bin, ensemblKnownSitesPath), out string vcfPath);
                proteinVariantDatabases.Add(
                    WriteSampleSpecificFasta(vcfPath, ensemblGenome, geneModelGtfOrGff, 7, Path.Combine(Path.GetDirectoryName(newBam), Path.GetFileNameWithoutExtension(newBam))));
            }
            STARWrapper.RemoveGenome(bin, genomeStarIndexDirectory);
        }

        public static string WriteSampleSpecificFasta(string vcfPath, Genome genome, string geneModelGtfOrGff, int minPeptideLength, string outprefix)
        {
            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).ToList();
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
            geneModel.AmendTranscripts(singleNucleotideVariants);
            List<Protein> proteins = new List<Protein>();
            List<Transcript> transcripts = geneModel.Genes.SelectMany(g => g.Transcripts).ToList();
            for (int i = 0; i < transcripts.Count; i++)
            {
                Console.WriteLine("Processing transcript " + i.ToString() + "/" + transcripts.Count.ToString() + " " + transcripts[i].ID + " " + transcripts[i].ProteinID);
                proteins.AddRange(transcripts[i].Translate(true, true));
            }
            int transcriptsWithVariants = geneModel.Genes.Sum(g => g.Transcripts.Count(x => x.CodingDomainSequences.Any(y => y.Variants.Count > 0)));
            string proteinVariantDatabase =  outprefix + ".protein.fasta";
            List<Protein> proteinsToWrite = proteins.OrderBy(p => p.Accession).Where(p => p.BaseSequence.Length >= minPeptideLength).ToList();
            ProteinDbWriter.WriteFastaDatabase(proteinsToWrite, proteinVariantDatabase, " ");
            using (StreamWriter writer = new StreamWriter(outprefix + ".protein.metrics"))
            {
                writer.WriteLine(Transcript.combinatoricFailures.ToString() + "\ttranscripts had too many heterozygous variants for combinatorics");
                Transcript.combinatoricFailures = 0;
                writer.WriteLine(proteinsToWrite.Count.ToString() + "\tprotein sequences");
                writer.WriteLine(proteinsToWrite.Min(p => p.BaseSequence.Length).ToString() + "\tminimum length");
                writer.WriteLine(proteinsToWrite.Max(p => p.BaseSequence.Length).ToString() + "\tmaxium length");
                writer.WriteLine(proteinsToWrite.Average(p => p.BaseSequence.Length).ToString() + "\taverage length");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf(ProteinAnnotation.SingleAminoAcidVariantLabel) > 0).ToString() + "\tSAV sequences");
                List<int> instancesOfSavs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace(ProteinAnnotation.SingleAminoAcidVariantLabel, "").Length) / ProteinAnnotation.SingleAminoAcidVariantLabel.Length).ToList();
                writer.WriteLine(instancesOfSavs.Max().ToString() + "\tmaximum SAVs in a sequence");
                writer.WriteLine(instancesOfSavs.Where(x => x > 0).Average().ToString() + "\taverage SAVs in sequence with one");
            }
            return proteinVariantDatabase;
        }
    }
}
