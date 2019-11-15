using Bio;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Gene :
        Interval
    {
        /// <summary>
        /// Constructing a gene using Bio object
        /// </summary>
        /// <param name="id"></param>
        /// <param name="chromosome"></param>
        /// <param name="strand"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="metadata"></param>
        public Gene(string id, Chromosome chromosome, string source, string strand, long oneBasedStart, long oneBasedEnd, MetadataListItem<List<string>> featureMetadata)
            : base(chromosome, chromosome.Sequence.ID, source, strand, oneBasedStart, oneBasedEnd)
        {
            ID = id;
            Chromosome = chromosome;
            FeatureMetadata = featureMetadata;
        }

        /// <summary>
        /// ID for this gene
        /// </summary>
        public string ID { get; set; }

        /// <summary>
        /// Chromosome where this gene is found (easier than accessing and casting parent every time)
        /// </summary>
        public Chromosome Chromosome { get; set; }

        /// <summary>
        /// Transcripts within this genez
        /// </summary>
        public List<Transcript> Transcripts { get; set; } = new List<Transcript>();

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "gene";

        public MetadataListItem<List<string>> FeatureMetadata { get; private set; }

        /// <summary>
        /// Translates all transcripts from this gene into protein sequences
        /// </summary>
        /// <param name="translateCodingDomains"></param>
        /// <param name="incompleteTranscriptAccessions"></param>
        /// <param name="selenocysteineContaining"></param>
        /// <returns></returns>
        public List<Protein> Translate(bool translateCodingDomains, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            incompleteTranscriptAccessions = incompleteTranscriptAccessions ?? new HashSet<string>();
            List<Protein> proteins = new List<Protein>();
            foreach (Transcript t in translateCodingDomains ? Transcripts.Where(t => t.IsProteinCoding()) : Transcripts)
            {
                if (incompleteTranscriptAccessions.Contains(t.ProteinID)) { continue; }
                lock (proteins) { proteins.Add(t.Protein(selenocysteineContaining)); }
            }
            return proteins;
        }

        public List<MetadataListItem<List<string>>> GetFeatures()
        {
            var features = new List<MetadataListItem<List<string>>>();
            var geneMetadata = GetGtfFeatureMetadata();
            features.Add(geneMetadata);
            features.AddRange(Transcripts.OrderBy(t => t.OneBasedStart).SelectMany(t => t.GetFeatures()));
            return features;
        }

        public override string GetGtfAttributes()
        {
            var attributes = GeneModel.SplitAttributes(FeatureMetadata.FreeText);
            List<Tuple<string, string>> attributeSubsections = new List<Tuple<string, string>>();

            string geneIdLabel = "gene_id";
            bool hasGeneId = attributes.TryGetValue(geneIdLabel, out string geneId);
            if (hasGeneId) { attributeSubsections.Add(new Tuple<string, string>(geneIdLabel, geneId)); }

            string geneNameLabel = "gene_name";
            bool hasGeneName = attributes.TryGetValue(geneNameLabel, out string geneName);
            if (hasGeneName) { attributeSubsections.Add(new Tuple<string, string>(geneNameLabel, geneName)); }

            string geneVersionLabel = "gene_version";
            bool hasGeneVersion = attributes.TryGetValue(geneVersionLabel, out string geneVersion);
            if (hasGeneVersion) { attributeSubsections.Add(new Tuple<string, string>(geneVersionLabel, geneVersion)); }

            string geneBiotypeLabel = "gene_biotype";
            bool hasGeneBiotype = attributes.TryGetValue(geneBiotypeLabel, out string geneBiotype);
            if (hasGeneBiotype) { attributeSubsections.Add(new Tuple<string, string>(geneBiotypeLabel, geneBiotype)); }

            return String.Join(" ", attributeSubsections.Select(x => x.Item1 + " \"" + x.Item2 + "\";"));
        }
    }
}