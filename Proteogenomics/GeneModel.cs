using Bio;
using Bio.IO.Gff;
using Bio.VCF;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace Proteogenomics
{
    /// <summary>
    /// Contains representation of genes, transcripts, exons, etc. represented in a gene model. Can be amended with variants.
    /// </summary>
    public class GeneModel
    {

        #region Private Fields

        /// <summary>
        /// Gets the first instance of a word
        /// </summary>
        private static Regex attributeKey = new Regex(@"([\w]+)");

        /// <summary>
        /// Gets anything inside quotes
        /// </summary>
        private static Regex attributeValue = new Regex(@"""([\w.]+)""");

        #endregion Private Fields

        #region Public Properties

        /// <summary>
        /// Genome this gene model is based on.
        /// </summary>
        public Genome Genome { get; set; }

        /// <summary>
        /// Genes represented in this gene model.
        /// </summary>
        public List<Gene> Genes { get; set; } = new List<Gene>();

        /// <summary>
        /// Start codons represented in this gene model.
        /// </summary>
        public List<Exon> StartCDS { get; set; } = new List<Exon>();

        #endregion Public Properties

        #region Public Constructor

        /// <summary>
        /// Constructs this GeneModel object from a Genome object and a corresponding GTF or GFF3 gene model file.
        /// </summary>
        /// <param name="genome"></param>
        /// <param name="geneModelFile"></param>
        public GeneModel(Genome genome, string geneModelFile)
        {
            this.Genome = genome;
            ReadGeneFeatures(geneModelFile);
        }

        #endregion Public Constructor

        #region Public Methods

        private Gene currentGene = null;
        private Transcript currentTranscript = null;

        public void ReadGeneFeatures(string geneModelFile)
        {
            ForceGffVersionTo2(geneModelFile, out string geneModelWithVersion2MarkedPath);
            List<ISequence> geneFeatures = new GffParser().Parse(geneModelWithVersion2MarkedPath).ToList();

            foreach (ISequence chromFeatures in geneFeatures)
            {
                ISequence chromSeq = Genome.Chromosomes.FirstOrDefault(x => x.ID.Split(' ')[0] == chromFeatures.ID);
                if (chromSeq == null) continue;

                chromFeatures.Metadata.TryGetValue("features", out object f);
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
                            key = attributeKey.Match(attrib.TrimStart()).Groups[1].Value;
                            val = attributeValue.Match(attrib.TrimStart()).Groups[1].Value;
                        }
                        if (!attributes.TryGetValue(key, out string x)) // sometimes there are two tags, so avoid adding twice
                            attributes.Add(key, val);
                    }

                    if (feature.FreeText.Contains('='))
                        ProcessGff3Feature(feature, start, end, chromSeq, attributes);
                    else
                        ProcessGtfFeature(feature, start, end, chromSeq, attributes);

                }
            }
        }

        /// <summary>
        /// Processes a feature from a GFF3 gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="OneBasedStart"></param>
        /// <param name="OneBasedEnd"></param>
        /// <param name="chromSeq"></param>
        /// <param name="attributes"></param>
        public void ProcessGff3Feature(MetadataListItem<List<string>> feature, long OneBasedStart, long OneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string geneId);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcriptId);
            bool hasTranscriptVersion = attributes.TryGetValue("version", out string transcriptVersion) && hasTranscriptId;
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string proteinId);

            if (hasGeneId && (currentGene == null || hasGeneId && geneId != currentGene.ID))
            {
                currentGene = new Gene(geneId, chromSeq, feature);
                Genes.Add(currentGene);
            }

            if (hasTranscriptId && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, feature);
                currentGene.Transcripts.Add(currentTranscript);
            }

            if (hasExonId || hasProteinId)
            {
                ISequence exon_dna = chromSeq.GetSubSequence(OneBasedStart - 1, OneBasedEnd - OneBasedStart + 1);

                Exon exon = new Exon(exon_dna, OneBasedStart, OneBasedEnd, chromSeq.ID, feature);

                if (hasExonId)
                {
                    //currentGene.Exons.Add(exon);
                    currentTranscript.Exons.Add(exon);
                }
                else if (hasProteinId)
                {
                    if (currentTranscript.CodingDomainSequences.Count == 0) StartCDS.Add(exon);
                    //currentGene.CodingDomainSequences.Add(exon);
                    currentTranscript.CodingDomainSequences.Add(exon);
                    currentTranscript.ProteinID = proteinId;
                }
            }
        }

        /// <summary>
        /// Processes a feature from a GTF gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="OneBasedStart"></param>
        /// <param name="OneBasedEnd"></param>
        /// <param name="chromSeq"></param>
        /// <param name="attributes"></param>
        public void ProcessGtfFeature(MetadataListItem<List<string>> feature, long OneBasedStart, long OneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string geneId);
            bool hasGeneName = attributes.TryGetValue("gene_name", out string geneName);
            bool hasGeneVersion = attributes.TryGetValue("gene_version", out string geneVersion);
            bool hasGeneBiotype = attributes.TryGetValue("gene_biotype", out string geneBiotype);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcriptId);
            bool hasTranscriptVersion = attributes.TryGetValue("transcript_version", out string transcriptVersion);
            bool hasTranscriptBiotype = attributes.TryGetValue("transcript_biotype", out string transcriptBiotype);
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasExonVersion = attributes.TryGetValue("exon_version", out string exonVersion);
            bool hasExonNumber = attributes.TryGetValue("exon_number", out string exonNumber);
            bool hasNearestRef = attributes.TryGetValue("nearest_ref", out string nearestRef); // Cufflinks
            bool hasClassCode = attributes.TryGetValue("class_code", out string classCode); // Cufflinks

            // Catch the transcript features before they go by if available, i.e. if the file doesn't just have exons
            if (feature.Key == "transcript" && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chromSeq, feature);
                    Genes.Add(currentGene);
                }

                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, feature);
                currentGene.Transcripts.Add(currentTranscript);
            }

            if (feature.Key == "exon" || feature.Key == "CDS")
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chromSeq, feature);
                    Genes.Add(currentGene);
                }

                if (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID)
                {
                    currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, null);
                    currentGene.Transcripts.Add(currentTranscript);
                }

                ISequence exon_dna = chromSeq.GetSubSequence(OneBasedStart - 1, OneBasedEnd - OneBasedStart + 1);

                Exon exon = new Exon(exon_dna, OneBasedStart, OneBasedEnd, chromSeq.ID, feature);

                if (feature.Key == "exon")
                {
                    //currentGene.Exons.Add(exon);
                    currentTranscript.Exons.Add(exon);
                }
                else if (feature.Key == "CDS")
                {
                    //currentGene.CodingDomainSequences.Add(exon);
                    currentTranscript.CodingDomainSequences.Add(exon);
                }
            }
        }

        /// <summary>
        /// Enters variant information into transcripts based on SnpEff annotations and into exons based on location.
        /// </summary>
        /// <param name="superVariants"></param>
        public void AmendTranscripts(List<VariantSuperContext> superVariants, string reference)
        {
            int binSize = 100000;
            Dictionary<Tuple<string, long>, List<VariantSuperContext>> chrIndexVariants = new Dictionary<Tuple<string, long>, List<VariantSuperContext>>();
            Dictionary<string, List<SnpEffAnnotation>> transcriptIdSnpEffVariants = new Dictionary<string, List<SnpEffAnnotation>>();
            foreach (VariantSuperContext superVariant in superVariants)
            {
                var key = new Tuple<string, long>(superVariant.Chr, superVariant.Start / binSize * binSize);
                if (chrIndexVariants.TryGetValue(key, out List<VariantSuperContext> vars))
                    vars.Add(superVariant);
                else
                    chrIndexVariants[key] = new List<VariantSuperContext> { superVariant };

                foreach (SnpEffAnnotation a in superVariant.SnpEffAnnotations)
                {
                    if (transcriptIdSnpEffVariants.TryGetValue(a.FeatureID, out List<SnpEffAnnotation> asdf))
                        asdf.Add(a);
                    else
                        transcriptIdSnpEffVariants.Add(a.FeatureID, new List<SnpEffAnnotation> { a });
                }
            }

            foreach (Transcript t in Genes.SelectMany(g => g.Transcripts))
            {
                string id = reference.StartsWith("GRCh38") ?
                    t.ID + "." + t.Version :
                    t.ID;
                if (transcriptIdSnpEffVariants.TryGetValue(id, out List<SnpEffAnnotation> annotations))
                    t.SnpEffVariants = annotations;
            }

            foreach (Exon x in Genes.SelectMany(g => g.Transcripts.SelectMany(t => t.Exons.Concat(t.CodingDomainSequences))))
            {
                for (long i = x.OneBasedStart / binSize; i < x.OneBasedEnd / binSize + 1; i++)
                {
                    var key = new Tuple<string, long>(x.ChromID.Split(' ')[0], i * binSize);
                    if (chrIndexVariants.TryGetValue(key, out List<VariantSuperContext> nearby_variants))
                        x.Variants = nearby_variants.Where(v => x.Includes(v.Start)).OfType<VariantContext>().ToList();
                }
            }
        }

        #endregion Public Methods

        #region Translation Methods

        public List<Protein> Translate(bool translateCodingDomains, bool translateWithVariants, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            return Genes.SelectMany(g => g.Translate(translateCodingDomains, translateWithVariants, incompleteTranscriptAccessions, selenocysteineContaining)).ToList();
        }

        public List<Protein> TranslateUsingAnnotatedStartCodons(GeneModel genesWithCodingDomainSequences, bool translateWithVariants, int minPeptideLength = 7)
        {
            return Genes.SelectMany(g => g.TranslateUsingAnnotatedStartCodons(genesWithCodingDomainSequences, translateWithVariants, minPeptideLength)).ToList();
        }

        #endregion Translation Methods

        #region Private Method

        private static Regex gffVersion = new Regex(@"(##gff-version\s+)(\d)");

        /// <summary>
        /// Required for using DotNetBio because it only handles GFF version 2 in the header.
        /// The only difference in the new version is within the attributes, which are stored as free text anyway.
        /// </summary>
        /// <param name="gffPath"></param>
        /// <param name="gffWithVersionMarked2Path"></param>
        private static void ForceGffVersionTo2(string gffPath, out string gffWithVersionMarked2Path)
        {
            gffWithVersionMarked2Path = Path.Combine(Path.GetDirectoryName(gffPath), Path.GetFileNameWithoutExtension(gffPath) + ".gff2" + Path.GetExtension(gffPath));
            if (File.Exists(gffWithVersionMarked2Path)) return;

            using (StreamReader reader = new StreamReader(gffPath))
            using (StreamWriter writer = new StreamWriter(gffWithVersionMarked2Path))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null)
                    {
                        break;
                    }
                    if (line.StartsWith("##gff-version"))
                    {
                        writer.Write(gffVersion.Replace(line, m => m.Groups[1] + "2") + "\n");
                    }
                    else
                    {
                        writer.Write(line + '\n');
                    }
                }
            }
        }

        #endregion Private Method

    }
}
