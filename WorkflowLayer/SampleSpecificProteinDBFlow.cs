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
    public class SampleSpecificProteinDBFlow
        : SpritzFlow
    {
        public const string Command = "proteins";
        public static readonly int MinimumPeptideLength = 7;

        public SampleSpecificProteinDBFlow()
            : base(MyWorkflow.SampleSpecificProteinDB)
        {
        }

        public SampleSpecificProteinDBParameters Parameters { get; set; }
        public List<string> ProteinVariantDatabases { get; private set; } = new List<string>();

        /// <summary>
        /// Generate sample specific protein database starting with fastq files
        /// </summary>
        public void GenerateSAVProteinsFromFastqs()
        {
            EnsemblDownloadsWrapper downloads = new EnsemblDownloadsWrapper();
            downloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            STARAlignmentFlow alignment = new STARAlignmentFlow();
            alignment.Parameters = new STARAlignmentParameters(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Reference, Parameters.Threads,
                Parameters.Fastqs, Parameters.StrandSpecific, Parameters.InferStrandSpecificity,
                Parameters.OverwriteStarAlignment, Parameters.GenomeStarIndexDirectory, downloads.ReorderedFastaPath,
                Parameters.ProteinFasta, Parameters.GeneModelGtfOrGff, Parameters.UseReadSubset, Parameters.ReadSubset);
            alignment.PerformTwoPassAlignment();
            downloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFasta);
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.SpritzDirectory, Parameters.GeneModelGtfOrGff, downloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(Parameters.SpritzDirectory, filteredGeneModelForScalpel);

            // Variant Calling
            VariantCallingFlow variantCalling = new VariantCallingFlow();
            variantCalling.CallVariants(Parameters.SpritzDirectory, Parameters.Reference, Parameters.Threads, Parameters.GenomeFasta, sortedBed12Path, Parameters.EnsemblKnownSitesPath, alignment.DedupedBamFiles, downloads.ReorderedFastaPath);

            // Generate databases
            GeneModel geneModel = new GeneModel(downloads.EnsemblGenome, Parameters.GeneModelGtfOrGff);
            ProteinVariantDatabases = variantCalling.AnnotatedVcfFilePaths.Select(annotatedVcf =>
                WriteSampleSpecificFasta(annotatedVcf, downloads.EnsemblGenome, geneModel, Parameters.Reference, downloads.ProteinAccessionSequence, downloads.BadProteinAccessions, downloads.SelenocysteineProteinAccessions, MinimumPeptideLength, Path.Combine(Path.GetDirectoryName(annotatedVcf), Path.GetFileNameWithoutExtension(annotatedVcf))))
                .ToList();
        }

        public string WriteSampleSpecificFasta(string vcfPath, Genome genome, GeneModel geneModel, string reference, Dictionary<string, string> proteinSeqeunces,
            HashSet<string> badProteinAccessions, Dictionary<string, string> selenocysteineContaininAccessions, int minPeptideLength, string outprefix)
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

        /// <summary>
        /// Writes out metrics regarding the combinitorial variants produced and entered into the fasta
        /// </summary>
        /// <param name="outprefix"></param>
        /// <param name="proteinsToWrite"></param>
        /// <returns></returns>
        public string WriteProteinFastaMetrics(string outprefix, List<Protein> proteinsToWrite)
        {
            string proteinMetrics = outprefix + ".variantprotein.fasta.metrics";
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

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (SampleSpecificProteinDBParameters)parameters;
            GenerateSAVProteinsFromFastqs();
        }
    }
}