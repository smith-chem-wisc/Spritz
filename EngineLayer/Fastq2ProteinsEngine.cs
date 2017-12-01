using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class Fastq2ProteinsEngine
    {
        public static void RunFromSra(string bin, string analysisDirectory, string reference, int threads, string sraAccession, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ucscKnownSitesPath, out string proteinVariantDatabase, bool useReadSubset = false, int readSubset = (int)10e6)
        {
            SRAToolkitWrapper.Fetch(bin, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
            RunFromFastqs(bin, analysisDirectory, reference, threads, fastqPaths, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, geneModelGtfOrGff, ucscKnownSitesPath, out proteinVariantDatabase, useReadSubset, readSubset);
        }

        public static void RunFromFastqs(string bin, string analysisDirectory, string reference, int threads, string[] fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, string ucscKnownSitesPath, out string proteinVariantDatabase, bool useReadSubset = false, int readSubset = (int)10e6)
        {
            if (Path.GetExtension(genomeFasta) == ".gz")
            {
                WrapperUtility.RunBashCommand("gunzip", WrapperUtility.ConvertWindowsPath(genomeFasta));
                genomeFasta = Path.ChangeExtension(genomeFasta, null);
            }

            SkewerWrapper.Trim(bin, threads, 19, fastqs, out string[] trimmedFastqs, out string skewerLog);

            Directory.CreateDirectory(genomeStarIndexDirectory);
            if (!File.Exists(Path.Combine(genomeStarIndexDirectory, "SA")))
            {
                STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { genomeFasta }, geneModelGtfOrGff);
            }

            string outPrefix = Path.Combine(Path.GetDirectoryName(trimmedFastqs[0]), Path.GetFileNameWithoutExtension(trimmedFastqs[0]));
            if (!File.Exists(outPrefix + STARWrapper.BamFileSuffix) || overwriteStarAlignment)
            {
                STARWrapper.LoadGenome(bin, genomeStarIndexDirectory);
                if (inferStrandSpecificity)
                {
                    STARWrapper.SubsetFastqs(trimmedFastqs, readSubset, bin, out string[] subsetFastqs);
                    if (useReadSubset) trimmedFastqs = subsetFastqs;
                    string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                    STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep);
                    strandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                }
                STARWrapper.BasicAlignReads(bin, threads, genomeStarIndexDirectory, trimmedFastqs, outPrefix, strandSpecific, STARGenomeLoadOption.LoadAndRemove);
            }

            GATKWrapper.PrepareBamAndFasta(bin, threads, outPrefix + STARWrapper.BamFileSuffix, genomeFasta, reference, out string newBam, out string ucscGenomeFasta);
            GATKWrapper.RealignIndels(bin, threads, ucscGenomeFasta, newBam, out string realignedBam); 
            GATKWrapper.VariantCalling(bin, threads, ucscGenomeFasta, realignedBam, Path.Combine(bin, ucscKnownSitesPath), out string vcfPath);

            VCFParser vcf = new VCFParser(vcfPath);
            List<VariantContext> singleNucleotideVariants = vcf.Select(x => x).Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1)).ToList();
            Genome genome = new Genome(genomeFasta);
            GeneModel geneModel = new GeneModel(genome, geneModelGtfOrGff);
            geneModel.AmendTranscripts(singleNucleotideVariants);
            List<Protein> proteins = geneModel.Genes.SelectMany(g => g.Transcripts.SelectMany(t => t.Translate(true, true))).ToList();
            proteinVariantDatabase = Path.Combine(Path.GetDirectoryName(fastqs[0]), Path.GetFileNameWithoutExtension(fastqs[0]) + ".protein.fasta");
            ProteinDbWriter.WriteFastaDatabase(proteins, proteinVariantDatabase, " ");
        }
    }
}
