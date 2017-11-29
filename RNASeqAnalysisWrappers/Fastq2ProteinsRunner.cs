using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace RNASeqAnalysisWrappers
{
    public class Fastq2ProteinsRunner
    {
        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, out string proteinVariantDatabase)
        {
            SRAToolkitWrapper.Install(bin);
            SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqPaths, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, geneModelGtfOrGff, out proteinVariantDatabase);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, string[] fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, out string proteinVariantDatabase)
        {
            SkewerWrapper.Install(bin);
            STARWrapper.Install(bin);
            RSeQCWrapper.Install(bin);
            GATKWrapper.Install(bin);

            if (Path.GetExtension(genomeFasta) == ".gz")
            {
                WrapperUtility.RunBashCommand("gunzip", WrapperUtility.ConvertWindowsPath(genomeFasta));
                genomeFasta = Path.ChangeExtension(genomeFasta, null);
            }

            SkewerWrapper.Trim(bin, Environment.ProcessorCount, 19, fastqs, out string[] trimmedFastqs, out string skewerLog);

            Directory.CreateDirectory(genomeStarIndexDirectory);
            if (!File.Exists(Path.Combine(genomeStarIndexDirectory, "SA")))
            {
                STARWrapper.GenerateGenomeIndex(bin, Environment.ProcessorCount, genomeStarIndexDirectory, new string[] { genomeFasta }, geneModelGtfOrGff);
            }

            string outPrefix = Path.Combine(Path.GetDirectoryName(trimmedFastqs[0]), Path.GetFileNameWithoutExtension(trimmedFastqs[0]));
            if (!File.Exists(outPrefix + STARWrapper.BamFileSuffix) || overwriteStarAlignment)
            {
                STARWrapper.LoadGenome(bin, genomeStarIndexDirectory);
                if (inferStrandSpecificity)
                {
                    STARWrapper.SubsetFastqs(trimmedFastqs, (int)10e6, bin, out string[] subsetFastqs);
                    string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                    STARWrapper.BasicAlignReads(bin, Environment.ProcessorCount, genomeStarIndexDirectory, trimmedFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep);
                    strandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff);
                }
                STARWrapper.BasicAlignReads(bin, Environment.ProcessorCount, genomeStarIndexDirectory, trimmedFastqs, outPrefix, strandSpecific, STARGenomeLoadOption.LoadAndRemove);
            }

            GATKWrapper.DownloadAndSortKnownVariantSitesForEnsembl(bin, bin, true, reference, genomeFasta, out string sortedKnownSitesFilename);
            GATKWrapper.PrepareBam(bin, Environment.ProcessorCount, outPrefix + STARWrapper.BamFileSuffix, genomeFasta, out string newBam);
            GATKWrapper.RealignIndels(bin, Environment.ProcessorCount, genomeFasta, newBam, out string realignedBam); 
            GATKWrapper.VariantCalling(bin, Environment.ProcessorCount, genomeFasta, realignedBam, Path.Combine(bin, sortedKnownSitesFilename), out string vcfPath);

            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantContext> variants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length)).ToList();
            Genome genome = new Genome(genomeFasta);
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
            geneModel.AmendTranscripts(variants);
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Transcripts.SelectMany(t => t.Translate(true, true))).ToList();
            proteinVariantDatabase = Path.Combine(Path.GetDirectoryName(fastqs[0]), Path.GetFileNameWithoutExtension(fastqs[0]) + ".protein.fasta");
            ProteinDbWriter.WriteFastaDatabase(proteins, proteinVariantDatabase, " ");
        }
    }
}
