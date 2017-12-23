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
    public class GeneModel
    {

        #region Private Fields

        private static Regex attributeKey = new Regex(@"([\w]+)"); // first instance of a word
        private static Regex attributeValue = new Regex(@"""([\w.]+)"""); // anything inside the quotes

        #endregion Private Fields

        #region Public Properties

        public Genome Genome { get; set; }
        public List<Gene> Genes { get; set; } = new List<Gene>();
        public List<Exon> StartCDS { get; set; } = new List<Exon>();

        #endregion Public Properties

        #region Public Constructor

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

        public void ProcessGff3Feature(MetadataListItem<List<string>> feature, long OneBasedStart, long OneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string gene_id);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcript_id);
            bool hasExonId = attributes.TryGetValue("exon_id", out string exon_id);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string protein_id);

            if (hasGeneId && (currentGene == null || hasGeneId && gene_id != currentGene.ID))
            {
                currentGene = new Gene(gene_id, chromSeq, feature);
                Genes.Add(currentGene);
            }

            if (hasTranscriptId && (currentTranscript == null || hasTranscriptId && transcript_id != currentTranscript.ID))
            {
                currentTranscript = new Transcript(transcript_id, currentGene, feature);
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
                    currentTranscript.ProteinID = protein_id;
                }
            }
        }

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

                currentTranscript = new Transcript(transcriptId, currentGene, feature);
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
                    currentTranscript = new Transcript(transcriptId, currentGene, null);
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

        public void AmendTranscripts(List<VariantContext> variants)
        {
            int binSize = 100000;
            Dictionary<Tuple<string, long>, List<VariantContext>> chrIndexVariants = new Dictionary<Tuple<string, long>, List<VariantContext>>();
            foreach(VariantContext variant in variants)
            {
                var key = new Tuple<string, long>(variant.Chr, variant.Start / binSize * binSize);
                if (chrIndexVariants.TryGetValue(key, out List<VariantContext> vars))
                    vars.Add(variant);
                else chrIndexVariants[key] = new List<VariantContext> { variant };
            }

            foreach (Exon x in Genes.SelectMany(g => g.Transcripts.SelectMany(t => t.Exons.Concat(t.CodingDomainSequences))))
            {
                for (long i = x.OneBasedStart / binSize; i < x.OneBasedEnd / binSize + 1; i++)
                {
                    var key = new Tuple<string, long>(x.ChromID.Split(' ')[0], i * binSize);
                    if (chrIndexVariants.TryGetValue(key, out List<VariantContext> nearby_variants))
                        x.Variants = nearby_variants.Where(v => x.Includes(v.Start)).ToList();
                }
            }
        }

        #endregion Public Methods

        #region Translation Methods

        public List<Protein> Translate(bool translateCodingDomains, bool translateWithVariants, HashSet<string> badProteinAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            return Genes.SelectMany(g => g.Translate(translateCodingDomains, translateWithVariants, badProteinAccessions, selenocysteineContaining)).ToList();
        }

        public List<Protein> TranslateUsingAnnotatedStartCodons(GeneModel genesWithCodingDomainSequences, bool translateWithVariants, int minPeptideLength = 7)
        {
            return Genes.SelectMany(g => g.TranslateUsingAnnotatedStartCodons(genesWithCodingDomainSequences, translateWithVariants, minPeptideLength)).ToList();
        }

        #endregion Translation Methods

        #region Private Method

        private static Regex gffVersion = new Regex(@"(##gff-version\s+)(\d)");
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
