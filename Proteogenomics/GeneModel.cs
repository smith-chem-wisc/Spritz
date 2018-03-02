using Bio;
using Bio.IO.Gff;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

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
        /// Forest of intervals by chromosome name
        /// </summary>
        public IntervalForest GenomeForest { get; set; } = new IntervalForest();

        /// <summary>
        /// Genes represented
        /// </summary>
        public List<Gene> Genes { get; set; } = new List<Gene>();

        /// <summary>
        /// Intergenic regions
        /// </summary>
        public List<Intergenic> Intergenics { get; set; } = new List<Intergenic>();

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

        #region Public Methods -- Read Gene Model File

        private Gene currentGene = null;
        private Transcript currentTranscript = null;

        public void ReadGeneFeatures(string geneModelFile)
        {
            ForceGffVersionTo2(geneModelFile, out string geneModelWithVersion2MarkedPath);
            List<ISequence> geneFeatures = new GffParser().Parse(geneModelWithVersion2MarkedPath).ToList();

            foreach (ISequence chromFeatures in geneFeatures)
            {
                ISequence chromSeq = Genome.Chromosomes.FirstOrDefault(x => x.FriendlyName == chromFeatures.ID).Sequence;
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
            CreateUTRsAndIntergenicRegions();
            Parallel.ForEach(Genes, gene => gene.TranscriptTree.Build(gene.Transcripts));
            GenomeForest.Build();
        }

        /// <summary>
        /// Processes a feature from a GFF3 gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chromSeq"></param>
        /// <param name="attributes"></param>
        public void ProcessGff3Feature(MetadataListItem<List<string>> feature, long oneBasedStart, long oneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string geneId);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcriptId);
            bool hasTranscriptVersion = attributes.TryGetValue("version", out string transcriptVersion) && hasTranscriptId;
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string proteinId);
            string strand = feature.SubItems["strand"][0];

            if (hasGeneId && (currentGene == null || hasGeneId && geneId != currentGene.ID))
            {
                currentGene = new Gene(geneId, chromSeq, strand, oneBasedStart, oneBasedEnd, feature);
                Genes.Add(currentGene);
                GenomeForest.Add(currentGene);
            }

            if (hasTranscriptId && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (hasExonId || hasProteinId)
            {
                if (hasExonId)
                {
                    ISequence exon_dna = chromSeq.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                    Exon exon = new Exon(exon_dna, oneBasedStart, oneBasedEnd, chromSeq.ID, strand);
                    currentTranscript.Exons.Add(exon);
                }
                else if (hasProteinId)
                {
                    CDS cds = new CDS(chromSeq.ID, strand, oneBasedStart, oneBasedEnd);
                    currentTranscript.CodingDomainSequences.Add(cds);
                    currentTranscript.ProteinID = proteinId;
                }
            }
        }

        /// <summary>
        /// Processes a feature from a GTF gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chromSeq"></param>
        /// <param name="attributes"></param>
        public void ProcessGtfFeature(MetadataListItem<List<string>> feature, long oneBasedStart, long oneBasedEnd, ISequence chromSeq, Dictionary<string, string> attributes)
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
            string strand = feature.SubItems["strand"][0];

            // Catch the transcript features before they go by if available, i.e. if the file doesn't just have exons
            if (feature.Key == "transcript" && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chromSeq, strand, oneBasedStart, oneBasedEnd, feature);
                    Genes.Add(currentGene);
                    GenomeForest.Add(currentGene);
                }

                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (feature.Key == "exon" || feature.Key == "CDS")
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chromSeq, strand, oneBasedStart, oneBasedEnd, feature);
                    Genes.Add(currentGene);
                    GenomeForest.Add(currentGene);
                }

                if (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID)
                {
                    currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd);
                    currentGene.Transcripts.Add(currentTranscript);
                    GenomeForest.Add(currentTranscript);
                }

                if (feature.Key == "exon")
                {
                    ISequence exon_dna = chromSeq.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                    Exon exon = new Exon(exon_dna, oneBasedStart, oneBasedEnd, chromSeq.ID, strand);
                    currentTranscript.Exons.Add(exon);
                }
                else if (feature.Key == "CDS")
                {
                    CDS cds = new CDS(chromSeq.ID, strand, oneBasedStart, oneBasedEnd);
                    currentTranscript.CodingDomainSequences.Add(cds);
                }
            }
        }

        #endregion Public Methods -- Read Gene Model File

        #region Private Method -- Read Gene Model File

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

        #endregion Private Method -- Read Gene Model File

        #region Public Methods

        /// <summary>
        /// Creates UTRs for transcripts and
        /// </summary>
        public void CreateUTRsAndIntergenicRegions()
        {
            foreach (IntervalTree it in GenomeForest.Forest.Values)
            {
                Gene previousPositiveStrandGene = null;
                Gene previousNegativeStrandGene = null;
                foreach (Gene gene in it.Intervals.OfType<Gene>().OrderBy(g => g.OneBasedStart))
                {
                    Gene previous = gene.Strand == "+" ? previousPositiveStrandGene : previousNegativeStrandGene;
                    Intergenic intergenic = previous == null ? null : new Intergenic(gene.ChromosomeID, gene.Strand, (gene.Strand == "+" ? previousPositiveStrandGene : previousNegativeStrandGene).OneBasedEnd + 1, gene.OneBasedStart - 1);

                    if (gene.Strand == "+")
                    {
                        previousPositiveStrandGene = gene;
                    }
                    if (gene.Strand == "-")
                    {
                        previousNegativeStrandGene = gene;
                    }

                    if (intergenic != null && intergenic.Length() <= 0)
                    {
                        GenomeForest.Add(intergenic);
                    }

                    Chromosome chromSeq = Genome.Chromosomes.FirstOrDefault(x => x.ChromosomeID == gene.ChromosomeID);
                    foreach (Transcript t in gene.Transcripts)
                    {
                        GenomeForest.Add(t.CreateIntrons());
                        GenomeForest.Add(t.CreateUTRs());
                        GenomeForest.Add(t.CreateUpDown(chromSeq));
                    }
                }
            }
        }

        /// <summary>
        /// Apply a list of variants to the intervals within this gene model
        /// </summary>
        /// <param name="variants"></param>
        public void ApplyVariants(List<Variant> variants)
        {
            List<Transcript> resultingTranscripts = new List<Transcript>();
            foreach (Variant v in variants.OrderByDescending(v => v.OneBasedStart).ToList())
            {
                List<Interval> intervals = GenomeForest.Forest[v.ChromosomeID].Stab(v.OneBasedStart);
                foreach (Interval i in intervals)
                {
                    i.ApplyVariant(v);
                }
            }
            foreach (Transcript t in Genes.SelectMany(g => g.Transcripts))
            {
                List<Variant> transcriptVariants = t.Variants.OrderByDescending(v => v.OneBasedStart).ToList(); // reversed, so that the coordinates of each successive variant is not changed
                resultingTranscripts.AddRange(t.ApplyVariantsCombinitorially(transcriptVariants));
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
    }
}