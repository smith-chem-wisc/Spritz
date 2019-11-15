using Bio;
using Bio.IO.Gff;
using Fclp;
using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;

namespace GtfSharp
{
    internal class GtfSharp
    {
        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to GtfSharp!");
            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.ReferenceGenome)
                .As('f', "reference_genome_fasta")
                .Required()
                .WithDescription("Reference genome from ensembl (FASTA).");

            p.Setup(arg => arg.CustomGeneModel)
                .As('g', "custom_gene_model")
                .Required()
                .WithDescription("Custom gene model (GTF, GFF2, GFF3).");

            p.Setup(arg => arg.ReferenceGeneModel)
                .As('r', "reference_gene_model")
                .Required()
                .WithDescription("Reference gene model from Ensembl (GTF, GFF2, GFF3).");

            p.SetupHelp("h", "help")
                .Callback(text => Console.WriteLine(text));

            var result = p.Parse(args);

            FilterGtfEntriesWithoutStrand(p.Object.CustomGeneModel, p.Object.ReferenceGenome, p.Object.ReferenceGeneModel);
        }

        /// <summary>
        /// Filters GTF or GFF entries that lack strand information
        /// Can filter also by zero abundance stringtie estimates
        /// Add CDS at the end
        /// </summary>
        /// <param name="gtfPath"></param>
        /// <param name="gtfOutPath"></param>
        public static void FilterGtfEntriesWithoutStrand(string gtfPath, string referenceGenomePath, string referenceGeneModelPath, bool filterEntriesWithZeroAbundanceStringtieEstimates = false)
        {
            var chromFeatures = GeneModel.SimplerParse(gtfPath);
            string filteredGtfPath = Path.Combine(Path.GetDirectoryName(gtfPath), Path.GetFileNameWithoutExtension(gtfPath) + ".filtered.gtf");
            using (var file = File.Create(filteredGtfPath))
            {
                var formatter = new GffFormatter();
                foreach (var chromISeq in chromFeatures)
                {
                    List<MetadataListItem<List<string>>> filteredFeatures = new List<MetadataListItem<List<string>>>();
                    bool isMetadata = chromISeq.Metadata.TryGetValue("features", out object featuresObj);
                    if (isMetadata)
                    {
                        bool okayTranscript = false;
                        var features = featuresObj as List<MetadataListItem<List<string>>>;
                        foreach (var feature in features)
                        {
                            if (!feature.SubItems.TryGetValue("strand", out List<string> strandish)) { continue; }
                            var attributes = GeneModel.SplitAttributes(feature.FreeText);
                            if (feature.Key == "transcript")
                            {
                                bool okayFpkm = !filterEntriesWithZeroAbundanceStringtieEstimates ||
                                    attributes.TryGetValue("FPKM", out string fpkm) && double.TryParse(fpkm, out double fpkmValue) && fpkmValue > 0;
                                bool okayTpm = !filterEntriesWithZeroAbundanceStringtieEstimates ||
                                    attributes.TryGetValue("TPM", out string tpm) && double.TryParse(tpm, out double tpmValue) && tpmValue > 0;
                                okayTranscript = okayFpkm && okayTpm;
                            }
                            if (okayTranscript)
                            {
                                filteredFeatures.Add(feature);
                            }
                        }
                    }
                    chromISeq.Metadata["features"] = filteredFeatures;
                }
                formatter.Format(file, chromFeatures);
            }
            Genome ensemblGenome = new Genome(referenceGenomePath);
            GeneModel newGeneModel = new GeneModel(ensemblGenome, filteredGtfPath);
            GeneModel referenceGeneModel = new GeneModel(ensemblGenome, referenceGeneModelPath);
            newGeneModel.CreateCDSFromAnnotatedStartCodons(referenceGeneModel);
            string filteredGtfWithCdsPath = Path.Combine(Path.GetDirectoryName(filteredGtfPath), Path.GetFileNameWithoutExtension(filteredGtfPath) + ".withcds.gtf");
            newGeneModel.PrintToGTF(filteredGtfWithCdsPath);
        }

        public class ApplicationArguments
        {
            public string CustomGeneModel { get; set; }
            public string ReferenceGeneModel { get; set; }
            public string ReferenceGenome { get; set; }
        }
    }
}