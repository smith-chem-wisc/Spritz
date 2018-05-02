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
    /// <summary>
    /// Create a database of proteins with single amino acid variants (SAVs) using GATK and SnpEff.
    /// </summary>
    public class SAVProteinDBFlow
    {
        /// <summary>
        /// Generate sample specific database starting with SRA accession number
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="reference"></param>
        /// <param name="threads"></param>
        /// <param name="sraAccession"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="overwriteStarAlignment"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="proteinFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="ensemblKnownSitesPath"></param>
        /// <param name="proteinVariantDatabases"></param>
        /// <param name="useReadSubset"></param>
        /// <param name="readSubset"></param>
        public static void GenerateSAVProteinsFromSra(
            string binDirectory, string analysisDirectory, string reference, int threads, string sraAccession,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases,
            bool useReadSubset = false, int readSubset = 300000)
        {
            List<string[]> fastqs = new List<string[]>();
            string[] sras = sraAccession.Split(',');
            foreach (string sra in sras)
            {
                SRAToolkitWrapper.Fetch(binDirectory, sra, analysisDirectory, out string[] fastqPaths, out string logPath);
                fastqs.Add(fastqPaths);
            }
            GenerateSAVProteinsFromFastqs(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, genomeFasta, proteinFasta, geneModelGtfOrGff, ensemblKnownSitesPath, out proteinVariantDatabases, useReadSubset, readSubset);
        }

        /// <summary>
        /// Generate sample specific protein database starting with fastq files
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="reference"></param>
        /// <param name="threads"></param>
        /// <param name="fastqs"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="inferStrandSpecificity"></param>
        /// <param name="overwriteStarAlignment"></param>
        /// <param name="genomeStarIndexDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="proteinFasta"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="ensemblKnownSitesPath"></param>
        /// <param name="proteinVariantDatabases"></param>
        /// <param name="useReadSubset"></param>
        /// <param name="readSubset"></param>
        public static void GenerateSAVProteinsFromFastqs(
            string binDirectory, string analysisDirectory, string reference, int threads, List<string[]> fastqs,
            bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory,
            string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, out List<string> proteinVariantDatabases,
            bool useReadSubset = false, int readSubset = 300000)
        {
            EnsemblDownloadsWrapper.PrepareEnsemblGenomeFasta(genomeFasta, out Genome ensemblGenome, out string reorderedFasta);
            STARAlignmentFlow.PerformTwoPassAlignment(binDirectory, analysisDirectory, reference, threads, fastqs, strandSpecific, inferStrandSpecificity, overwriteStarAlignment, genomeStarIndexDirectory, reorderedFasta, proteinFasta, geneModelGtfOrGff, out List<string> firstPassSpliceJunctions, out string secondPassGenomeDirectory, out List<string> sortedBamFiles, out List<string> dedupedBamFiles, out List<string> chimericSamFiles, out List<string> chimericJunctionFiles, useReadSubset, readSubset);
            EnsemblDownloadsWrapper.GetImportantProteinAccessions(binDirectory, proteinFasta, out var proteinSequences, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineContainingAccessions);
            EnsemblDownloadsWrapper.FilterGeneModel(binDirectory, geneModelGtfOrGff, ensemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(binDirectory, filteredGeneModelForScalpel);

            // Variant Calling
            string scriptName = Path.Combine(binDirectory, "scripts", "variantCalling.bash");
            List<string> variantCallingCommands = new List<string>();
            List<string> vcfFilePaths = new List<string>();
            List<string> annotatedVcfFilePaths = new List<string>();
            List<string> snpEffHtmlFilePaths = new List<string>();
            List<string> annotatedGenesSummaryPaths = new List<string>();
            List<string> scapelVcfFilePaths = new List<string>();
            List<string> annotatedScapelVcfFilePaths = new List<string>();
            List<string> scalpelSnpEffHtmlFilePaths = new List<string>();
            List<string> scalpelAnnotatedGenesSummaryPaths = new List<string>();
            foreach (string dedupedBam in dedupedBamFiles)
            {
                // GATK
                variantCallingCommands.AddRange(GATKWrapper.SplitNCigarReads(binDirectory, genomeFasta, dedupedBam, out string splitTrimBam));
                variantCallingCommands.AddRange(GATKWrapper.VariantCalling(binDirectory, threads, reorderedFasta, splitTrimBam, Path.Combine(binDirectory, ensemblKnownSitesPath), out string vcfPath));
                vcfFilePaths.Add(vcfPath);
                variantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(binDirectory, reference, vcfPath, out string htmlReport, out string annotatedVcfPath, out string annotatedGenesSummaryPath));
                annotatedVcfFilePaths.Add(annotatedVcfPath);
                snpEffHtmlFilePaths.Add(htmlReport);
                annotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath);

                // Scalpel
                variantCallingCommands.AddRange(ScalpelWrapper.CallIndels(binDirectory, threads, genomeFasta, sortedBed12Path, splitTrimBam, Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + "_scalpelOut"), out string scalpelVcf));
                scapelVcfFilePaths.Add(scalpelVcf);
                variantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(binDirectory, reference, scalpelVcf, out string htmlReport2, out string annotatedScalpelVcfPath, out string annotatedGenesSummaryPath2));
                annotatedScapelVcfFilePaths.Add(annotatedScalpelVcfPath);
                scalpelSnpEffHtmlFilePaths.Add(htmlReport2);
                scalpelAnnotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath2);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, variantCallingCommands).WaitForExit();

            // Generate databases
            GeneModel geneModel = new GeneModel(ensemblGenome, geneModelGtfOrGff);
            proteinVariantDatabases = annotatedVcfFilePaths.Select(annotatedVcf =>
                WriteSampleSpecificFasta(annotatedVcf, ensemblGenome, geneModel, reference, proteinSequences, badProteinAccessions, selenocysteineContainingAccessions, 7, Path.Combine(Path.GetDirectoryName(annotatedVcf), Path.GetFileNameWithoutExtension(annotatedVcf))))
                .ToList();
        }

        public static string WriteSampleSpecificFasta(string vcfPath, Genome genome, GeneModel geneModel, string reference, Dictionary<string, string> proteinSeqeunces, HashSet<string> badProteinAccessions, Dictionary<string, string> selenocysteineContaininAccessions, int minPeptideLength, string outprefix)
        {
            if (!File.Exists(vcfPath))
            {
                Console.WriteLine("Error: VCF not found: " + vcfPath);
                return "Error: VCF not found: " + vcfPath;
            }

            // Parse VCF file
            VCFParser vcf = new VCFParser(vcfPath);
            List<Variant> singleNucleotideVariants = vcf.Select(x => x)
                .Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1))
                .Select(v => new Variant(null, v, genome)).ToList();

            // Apply the variants combinitorially, and translate the variant transcripts
            List<Transcript> variantTranscripts = geneModel.ApplyVariants(singleNucleotideVariants);
            List<Protein> variantProteins = new List<Protein>();
            for (int i = 0; i < variantTranscripts.Count; i++)
            {
                if (badProteinAccessions.Contains(variantTranscripts[i].ID) || badProteinAccessions.Contains(variantTranscripts[i].ProteinID))
                {
                    continue;
                }
                variantProteins.Add(variantTranscripts[i].Protein(selenocysteineContaininAccessions));
            }

            int transcriptsWithVariants = variantTranscripts.Count(t => t.CodingDomainSequences.Any(y => y.Variants.Count > 0));
            string proteinVariantDatabase = outprefix + ".variantprotein.fasta";
            string proteinVariantDatabaseXml = outprefix + ".variantprotein.xml";
            List<Protein> variantProteinsToWrite = variantProteins.OrderBy(p => p.Accession).Where(p => p.BaseSequence.Length >= minPeptideLength).ToList();
            ProteinDbWriter.WriteFastaDatabase(variantProteinsToWrite, proteinVariantDatabase, "|");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), variantProteinsToWrite, proteinVariantDatabaseXml);
            WriteProteinFastaMetrics(outprefix, variantProteinsToWrite);

            // Apply the SNVs as annotations, and translate the variant transcripts
            //List<Transcript> transcripts = geneModel.ApplyVariants(singleNucleotideVariants);
            //List<Protein> proteins = new List<Protein>();
            //for (int i = 0; i < transcripts.Count; i++)
            //{
            //    if (badProteinAccessions.Contains(transcripts[i].ID) || badProteinAccessions.Contains(transcripts[i].ProteinID))
            //    {
            //        continue;
            //    }
            //    proteins.Add(transcripts[i].Protein(selenocysteineContaininAccessions));
            //}
            //string proteinDatabaseXml = outprefix + ".protein.xml";
            //List<Protein> proteinsToWrite = proteins.OrderBy(p => p.Accession).Where(p => p.BaseSequence.Length >= minPeptideLength).ToList();
            //ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsToWrite, proteinDatabaseXml);

            return proteinVariantDatabase;
        }

        public static string WriteProteinFastaMetrics(string outprefix, List<Protein> proteinsToWrite)
        {
            string proteinMetrics = outprefix + ".protein.metrics";
            using (StreamWriter writer = new StreamWriter(proteinMetrics))
            {
                Transcript.combinatoricFailures = new List<string>();
                writer.WriteLine(proteinsToWrite.Count.ToString() + "\tprotein sequences");
                writer.WriteLine(proteinsToWrite.Min(p => p.BaseSequence.Length).ToString() + "\tminimum length");
                writer.WriteLine(proteinsToWrite.Max(p => p.BaseSequence.Length).ToString() + "\tmaxium length");
                writer.WriteLine(proteinsToWrite.Average(p => p.BaseSequence.Length).ToString() + "\taverage length");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf(FunctionalClass.MISSENSE.ToString()) > 0).ToString() + "\tSAV sequences");
                List<int> instancesOfSavs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace(FunctionalClass.MISSENSE.ToString(), "").Length) / FunctionalClass.MISSENSE.ToString().Length).ToList();
                writer.WriteLine(instancesOfSavs.Max().ToString() + "\tmaximum SAVs in a sequence");
                if (instancesOfSavs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSavs.Where(x => x > 0).Average().ToString() + "\taverage SAVs in sequence with one");
                writer.WriteLine();
                writer.WriteLine(proteinsToWrite.Count(p => p.FullName.IndexOf(FunctionalClass.SILENT.ToString()) > 0).ToString() + "\tsequences with synonymous codon variation");
                List<int> instancesOfSynonymousSnvs = proteinsToWrite.Select(p => (p.FullName.Length - p.FullName.Replace(FunctionalClass.SILENT.ToString(), "").Length) / FunctionalClass.SILENT.ToString().Length).ToList();
                writer.WriteLine(instancesOfSynonymousSnvs.Max().ToString() + "\tmaximum synonymous variations in a sequence");
                if (instancesOfSynonymousSnvs.Count(x => x > 0) > 0) writer.WriteLine(instancesOfSynonymousSnvs.Where(x => x > 0).Average().ToString() + "\taverage synonymous variations in sequence with one");

                writer.WriteLine();
                writer.WriteLine("Skipped due to combinitoric explosion (> 5 heterozygous nonsynonymous variations):");
                writer.WriteLine(Transcript.combinatoricFailures.Count.ToString() + "\theterozygous nonsynonymous variations skipped");
                foreach (string failure in Transcript.combinatoricFailures)
                {
                    writer.WriteLine(failure);
                }
            }
            return proteinMetrics;
        }    
    }
}