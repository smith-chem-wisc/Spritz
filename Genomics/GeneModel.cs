using System.Text.RegularExpressions;
using System.Linq;
using System.Collections.Generic;
using Bio;
using Bio.IO.Gff;
using Bio.VCF;

namespace Genomics
{
    public class GeneModel
    {

        #region Private Fields

        private static Regex attribute_key = new Regex(@"([\w]+)"); // first instance of a word
        private static Regex attribute_value = new Regex(@"""([\w.]+)"""); // anything inside the quotes

        #endregion Private Fields

        #region Public Properties

        public Genome genome { get; set; }
        public List<Gene> genes { get; set; } = new List<Gene>();
        public List<Exon> start_cds { get; set; } = new List<Exon>();

        #endregion Public Properties

        #region Public Constructor

        public GeneModel(Genome genome, string geneModelFile)
        {
            this.genome = genome;
            ReadGeneFeatures(geneModelFile);
        }

        #endregion Public Constructor

        #region Public Method

        private Gene current_gene = null;
        private Transcript current_transcript = null;

        public void ReadGeneFeatures(string geneModelFile)
        {
            List<ISequence> gene_features = new GffParser().Parse(geneModelFile).ToList();

            foreach (ISequence chrom_features in gene_features)
            {
                ISequence chromSeq = genome.chroms.FirstOrDefault(x => x.ID.Split(' ')[0] == chrom_features.ID);
                if (chromSeq == null) continue;

                chrom_features.Metadata.TryGetValue("features", out object f);
                List<MetadataListItem<List<string>>> features = f as List<MetadataListItem<List<string>>>;
                for (int i = 0; i < features.Count; i++)
                {
                    MetadataListItem<List<string>> feature = features[i];
                    long.TryParse(feature.SubItems["start"][0], out long start);
                    long.TryParse(feature.SubItems["end"][0], out long end);

                    Dictionary<string, string> attributes = new Dictionary<string, string>();
                    foreach (string attrib in feature.FreeText.Split(';'))
                    {
                        string key;
                        string val;
                        if (feature.FreeText.Contains('=')) // GFF3
                        {
                            key = attrib.Split('=')[0].TrimStart();
                            val = attrib.Split('=')[1].TrimStart();
                        }
                        else // GFF1 or GTF
                        {
                            key = attribute_key.Match(attrib.TrimStart()).Groups[1].Value;
                            val = attribute_value.Match(attrib.TrimStart()).Groups[1].Value;
                        }
                        if (!attributes.TryGetValue(key, out string x)) // sometimes there are two tags, so avoid adding twice
                            attributes.Add(key, val);
                    }

                    if (feature.FreeText.Contains('='))
                        process_gff3_feature(feature, start, end, chromSeq, attributes);
                    else
                        process_gtf_feature(feature, start, end, chromSeq, attributes);

                }
            }
        }

        public void process_gff3_feature(MetadataListItem<List<string>> feature, long OneBasedStart, long OneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool has_gene_id = attributes.TryGetValue("gene_id", out string gene_id);
            bool has_transcript_id = attributes.TryGetValue("transcript_id", out string transcript_id);
            bool has_exon_id = attributes.TryGetValue("exon_id", out string exon_id);
            bool has_protein_id = attributes.TryGetValue("protein_id", out string protein_id);

            if (has_gene_id && (current_gene == null || has_gene_id && gene_id != current_gene.ID))
            {
                current_gene = new Gene(gene_id, chromSeq.ID, feature);
                genes.Add(current_gene);
            }

            if (has_transcript_id && (current_transcript == null || has_transcript_id && transcript_id != current_transcript.ID))
            {
                current_transcript = new Transcript(transcript_id, current_gene, feature);
                current_gene.transcripts.Add(current_transcript);
            }

            if (has_exon_id || has_protein_id)
            {
                ISequence exon_dna = chromSeq.GetSubSequence(OneBasedStart - 1, OneBasedEnd - OneBasedStart + 1);

                Exon exon = new Exon(exon_dna, OneBasedStart, OneBasedEnd, chromSeq.ID, feature);

                if (has_exon_id)
                {
                    current_gene.exons.Add(exon);
                    current_transcript.Exons.Add(exon);
                }
                else if (has_protein_id)
                {
                    if (current_transcript.CDS.Count == 0) start_cds.Add(exon);
                    current_transcript.CDS.Add(exon);
                    current_transcript.ProteinID = protein_id;
                }
            }
        }

        public void process_gtf_feature(MetadataListItem<List<string>> feature, long OneBasedStart, long OneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool has_gene_id = attributes.TryGetValue("gene_id", out string gene_id);
            bool has_gene_name = attributes.TryGetValue("gene_name", out string gene_name);
            bool has_gene_version = attributes.TryGetValue("gene_version", out string gene_version);
            bool has_gene_biotype = attributes.TryGetValue("gene_biotype", out string gene_biotype);
            bool has_transcript_id = attributes.TryGetValue("transcript_id", out string transcript_id);
            bool has_transcript_version = attributes.TryGetValue("transcript_version", out string transcript_version);
            bool has_transcript_biotype = attributes.TryGetValue("transcript_biotype", out string transcript_biotype);
            bool has_exon_id = attributes.TryGetValue("exon_id", out string exon_id);
            bool has_exon_version = attributes.TryGetValue("exon_version", out string exon_version);
            bool has_exon_number = attributes.TryGetValue("exon_number", out string exon_number);
            bool has_nearest_ref = attributes.TryGetValue("nearest_ref", out string nearest_ref); // Cufflinks
            bool has_class_code = attributes.TryGetValue("class_code", out string class_code); // Cufflinks

            // Catch the transcript features before they go by if available, i.e. if the file doesn't just have exons
            if (feature.Key == "transcript" && (current_transcript == null || has_transcript_id && transcript_id != current_transcript.ID))
            {
                if (current_gene == null || has_gene_id && gene_id != current_gene.ID)
                {
                    current_gene = new Gene(gene_id, chromSeq.ID, feature);
                    genes.Add(current_gene);
                }

                current_transcript = new Transcript(transcript_id, current_gene, feature);
                current_gene.transcripts.Add(current_transcript);
            }

            if (feature.Key == "exon" || feature.Key == "CDS")
            {
                if (current_gene == null || has_gene_id && gene_id != current_gene.ID)
                {
                    current_gene = new Gene(gene_id, chromSeq.ID, feature);
                    genes.Add(current_gene);
                }

                if (current_transcript == null || has_transcript_id && transcript_id != current_transcript.ID)
                {
                    current_transcript = new Transcript(transcript_id, current_gene, null);
                    current_gene.transcripts.Add(current_transcript);
                }

                ISequence exon_dna = chromSeq.GetSubSequence(OneBasedStart - 1, OneBasedEnd - OneBasedStart + 1);

                Exon exon = new Exon(exon_dna, OneBasedStart, OneBasedEnd, chromSeq.ID, feature);

                if (feature.Key == "exon")
                {
                    current_gene.exons.Add(exon);
                    current_transcript.Exons.Add(exon);
                }
                else if (feature.Key == "CDS")
                {
                    current_transcript.CDS.Add(exon);
                }
            }
        }

        public void amend_transcripts(List<VariantContext> variants)
        {

        }

        #endregion Public Method

    }
}
