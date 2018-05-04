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

        public SampleSpecificProteinDBFlow()
            : base(MyWorkflow.SampleSpecificProteinDB)
        {
        }

        public SampleSpecificProteinDBParameters Parameters { get; set; }
        public List<string> IndelAppliedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> IndelAppliedProteinXmlDatabases { get; private set; } = new List<string>();
        public List<string> VariantAnnotatedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> VariantAnnotatedProteinXmlDatabases { get; private set; } = new List<string>();
        public List<string> VariantAppliedProteinFastaDatabases { get; private set; } = new List<string>();
        public List<string> VariantAppliedProteinXmlDatabases { get; private set; } = new List<string>();
        private EnsemblDownloadsWrapper Downloads { get; set; } = new EnsemblDownloadsWrapper();


        /// <summary>
        /// Generate sample specific protein database starting with fastq files
        /// </summary>
        public void GenerateSAVProteinsFromFastqs()
        {
            // Download references and align reads
            Downloads.PrepareEnsemblGenomeFasta(Parameters.GenomeFasta);
            STARAlignmentFlow alignment = new STARAlignmentFlow();
            alignment.Parameters = new STARAlignmentParameters(
                Parameters.SpritzDirectory, 
                Parameters.AnalysisDirectory, 
                Parameters.Reference, 
                Parameters.Threads,
                Parameters.Fastqs, 
                Parameters.StrandSpecific,
                Parameters.InferStrandSpecificity,
                Parameters.OverwriteStarAlignment, 
                Parameters.GenomeStarIndexDirectory, 
                Downloads.ReorderedFastaPath,
                Parameters.ReferenceGeneModelGtfOrGff,
                Parameters.UseReadSubset, 
                Parameters.ReadSubset);
            alignment.PerformTwoPassAlignment();
            Downloads.GetImportantProteinAccessions(Parameters.SpritzDirectory, Parameters.ProteinFasta);
            EnsemblDownloadsWrapper.FilterGeneModel(Parameters.SpritzDirectory, Parameters.ReferenceGeneModelGtfOrGff, Downloads.EnsemblGenome, out string filteredGeneModelForScalpel);
            string sortedBed12Path = BEDOPSWrapper.Gtf2Bed12(Parameters.SpritzDirectory, filteredGeneModelForScalpel);
            GeneModel referenceGeneModel = new GeneModel(Downloads.EnsemblGenome, Parameters.ReferenceGeneModelGtfOrGff);

            // Merge new and reference gene models, if a new one is specified
            Merge(referenceGeneModel, Parameters.NewGeneModelGtfOrGff);

            // Variant Calling
            VariantCallingFlow variantCalling = new VariantCallingFlow();
            variantCalling.CallVariants(
                Parameters.SpritzDirectory,
                Parameters.Reference, 
                Parameters.Threads, 
                Parameters.GenomeFasta,
                sortedBed12Path, 
                Parameters.EnsemblKnownSitesPath, 
                alignment.DedupedBamFiles, 
                Downloads.ReorderedFastaPath);

            // Apply Variations
            for (int i = 0; i < variantCalling.GatkVcfFilePaths.Count; i++)
            {
                // Indel database
                //string indelVcf = variantCalling.ScalpelAnnotatedVcfFilePaths[i];
                //var indelTranscripts = variantCalling.ApplyIndels(indelVcf, Downloads.EnsemblGenome, referenceGeneModel);
                //string indelOutPrefix = Path.Combine(Path.GetDirectoryName(indelVcf), Path.GetFileNameWithoutExtension(indelVcf) + "_indel");
                //(string, string, List<Protein>) indelDatabases = WriteSampleSpecificFasta(
                //    indelTranscripts,
                //    referenceGeneModel,
                //    indelOutPrefix);
                //IndelAppliedProteinFastaDatabases.Add(indelDatabases.Item1);
                //IndelAppliedProteinXmlDatabases.Add(indelDatabases.Item2);

                // Variant annotated database
                string snvVcf = variantCalling.GatkAnnotatedVcfFilePaths[i];
                variantCalling.AnnotateSAVs(snvVcf, Downloads.EnsemblGenome, referenceGeneModel);
                string annotatedOutPrefix = Path.Combine(Path.GetDirectoryName(snvVcf), Path.GetFileNameWithoutExtension(snvVcf) + "_snvAnnotated");
                (string, string, List<Protein>) annotatedDatabases = WriteSampleSpecificFasta(
                    referenceGeneModel.Genes.SelectMany(g => g.Transcripts).ToList(), 
                    referenceGeneModel,  
                    annotatedOutPrefix);
                VariantAnnotatedProteinFastaDatabases.Add(annotatedDatabases.Item1);
                VariantAnnotatedProteinXmlDatabases.Add(annotatedDatabases.Item2);

                // Variant applied database
                var variantTranscripts = variantCalling.ApplySNVs(snvVcf, Downloads.EnsemblGenome, referenceGeneModel);
                string appliedOutPrefix = Path.Combine(Path.GetDirectoryName(snvVcf), Path.GetFileNameWithoutExtension(snvVcf) + "_snvApplied");
                (string, string, List<Protein>) appliedDatabases = WriteSampleSpecificFasta(
                    variantTranscripts,
                    referenceGeneModel,
                    appliedOutPrefix);
                VariantAnnotatedProteinFastaDatabases.Add(appliedDatabases.Item1);
                VariantAnnotatedProteinXmlDatabases.Add(appliedDatabases.Item2);
                WriteProteinFastaMetrics(appliedOutPrefix, appliedDatabases.Item3);
            }
        }

        /// <summary>
        /// Read in a new gene model and merge it with this one
        /// </summary>
        /// <param name="alternativeGeneModelPath"></param>
        private void Merge(GeneModel referenceGeneModel, string alternativeGeneModelPath)
        {
            GeneModel newGeneModel = null;
            if (alternativeGeneModelPath != null && File.Exists(alternativeGeneModelPath))
            {
                string newGeneModelPath = EnsemblDownloadsWrapper.ConvertFirstColumnEnsembl2UCSC(Parameters.SpritzDirectory, Parameters.Reference, Parameters.NewGeneModelGtfOrGff);
                newGeneModel = new GeneModel(Downloads.EnsemblGenome, newGeneModelPath);
                newGeneModel.CreateCDSFromAnnotatedStartCodons(referenceGeneModel);
            }
            if (newGeneModel != null)
            {
                referenceGeneModel.Merge(newGeneModel);
            }
        }

        /// <summary>
        /// Write transcripts to protein fasta and xml databases
        /// </summary>
        /// <param name="transcripts"></param>
        /// <param name="genome"></param>
        /// <param name="geneModel"></param>
        /// <param name="badProteinAccessions"></param>
        /// <param name="selenocysteineContaininAccessions"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="outprefix"></param>
        /// <returns></returns>
        private (string, string, List<Protein>) WriteSampleSpecificFasta(List<Transcript> transcripts, GeneModel geneModel, string outprefix)
        {
            // Apply the variants combinitorially, and translate the variant transcripts
            List<Protein> variantProteins = transcripts
                .Where(t => !Downloads.BadProteinAccessions.Contains(t.ID) && !Downloads.BadProteinAccessions.Contains(t.ProteinID))
                .Select(t => t.Protein(Downloads.SelenocysteineProteinAccessions)).ToList();

            string proteinFasta = outprefix + ".variantprotein.fasta";
            string proteinXml = outprefix + ".variantprotein.xml";
            List<Protein> variantProteinsToWrite = variantProteins.OrderBy(p => p.Accession).Where(p => p.BaseSequence.Length >= Parameters.MinPeptideLength).ToList();
            ProteinDbWriter.WriteFastaDatabase(variantProteinsToWrite, proteinFasta, "|");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), variantProteinsToWrite, proteinXml);
            return (proteinFasta, proteinXml, variantProteins);
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