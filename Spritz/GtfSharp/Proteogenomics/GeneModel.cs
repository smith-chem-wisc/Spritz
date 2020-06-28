using Bio;
using Bio.Extensions;
using Bio.IO.Gff;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Proteogenomics
{
    /// <summary>
    /// Contains representation of genes, transcripts, exons, etc. represented in a gene model. Can be amended with variants.
    /// </summary>
    public class GeneModel
    {
        public const int CODON_SIZE = 3;
        private Gene currentGene = null;
        private Transcript currentTranscript = null;

        /// <summary>
        /// Gets the first instance of a word
        /// </summary>
        private static Regex AttributeKey = new Regex(@"([\w]+)");

        /// <summary>
        /// Gets anything inside quotes
        /// </summary>
        private static Regex AttributeValue = new Regex(@"""([^""]+)""");
        
        /// <summary>
        /// Constructs this GeneModel object from a Genome object and a corresponding GTF or GFF3 gene model file.
        /// </summary>
        /// <param name="genome"></param>
        /// <param name="geneModelFile"></param>
        public GeneModel(Genome genome, string geneModelFile)
        {
            Genome = genome;
            ReadGeneFeatures(geneModelFile);
        }

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

        #region Methods -- Read Gene Model File

        /// <summary>
        /// Reads gene model features into data structures contained within this library
        /// </summary>
        /// <param name="geneModelFile"></param>
        public void ReadGeneFeatures(string geneModelFile)
        {
            foreach (ISequence chromFeatures in SimplerParse(geneModelFile))
            {
                Chromosome chrom = Genome.Chromosomes.FirstOrDefault(x => x.FriendlyName == chromFeatures.ID);
                if (chrom == null) { continue; }
                chromFeatures.Metadata.TryGetValue("features", out object f);
                List<MetadataListItem<List<string>>> features = f as List<MetadataListItem<List<string>>>;
                for (int i = 0; i < features.Count; i++)
                {
                    MetadataListItem<List<string>> feature = features[i];
                    long.TryParse(feature.SubItems["start"][0], out long start);
                    long.TryParse(feature.SubItems["end"][0], out long end);
                    var attributes = SplitAttributes(feature.FreeText);
                    if (feature.FreeText.Contains('='))
                    {
                        ProcessGff3Feature(feature, start, end, chrom, attributes);
                    }
                    else
                    {
                        ProcessGtfFeature(feature, start, end, chrom, attributes);
                    }
                }
            }
            if (currentTranscript != null)
            {
                Transcript.SetRegions(currentTranscript);
                currentTranscript.FrameCorrection();
            }
            CreateIntergenicRegions();
            // possibly check transcript sanity here with Parallel.ForEach(Genes.SelectMany(g => g.Transcripts).ToList(), t => t.SanityCheck());
            GenomeForest.Build();
        }

        /// <summary>
        /// Processes a feature from a GFF3 gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chrom"></param>
        /// <param name="attributes"></param>
        public void ProcessGff3Feature(MetadataListItem<List<string>> feature, long oneBasedStart, long oneBasedEnd, Chromosome chrom, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string geneId);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcriptId);
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string proteinId);
            bool hasSource = feature.SubItems.TryGetValue("source", out List<string> sourceish); // false if empty ("." in GFF format)
            bool hasStrand = feature.SubItems.TryGetValue("strand", out List<string> strandish); // false if empty ("." in GFF format)
            bool hasFrame = feature.SubItems.TryGetValue("frame", out List<string> framey); // false if empty ("." in GFF format)

            string source = hasSource ? sourceish[0] : "";
            if (!hasStrand) { return; } // strand is a required to do anything in this program
            string strand = strandish[0];
            int frame = 0;
            if (hasFrame) { int.TryParse(framey[0], out frame); }

            if (hasGeneId && (currentGene == null || hasGeneId && geneId != currentGene.ID))
            {
                currentGene = new Gene(geneId, chrom, source, strand, oneBasedStart, oneBasedEnd, feature);
                Genes.Add(currentGene);
                GenomeForest.Add(currentGene);
            }

            if (hasTranscriptId && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                if (currentTranscript != null)
                {
                    Transcript.SetRegions(currentTranscript);
                    currentTranscript.FrameCorrection();
                }
                currentTranscript = new Transcript(transcriptId, currentGene, source, strand, oneBasedStart, oneBasedEnd, null, feature);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (hasExonId)
            {
                ISequence exon_dna = chrom.Sequence.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                Exon exon = new Exon(currentTranscript, currentTranscript.IsStrandPlus() ? exon_dna : exon_dna.GetReverseComplementedSequence(),
                    source, oneBasedStart, oneBasedEnd, chrom == null ? "" : chrom.ChromosomeID, strand, feature);
                if (exon.Length() > 0)
                {
                    currentTranscript.Exons.Add(exon);
                }
            }
            else if (hasProteinId)
            {
                CDS cds = new CDS(currentTranscript, chrom.Sequence.ID, source, strand, oneBasedStart, oneBasedEnd, frame);
                if (cds.Length() > 0)
                {
                    currentTranscript.CodingDomainSequences.Add(cds);
                    currentTranscript.ProteinID = proteinId;
                }
            }
            else // nothing to do
            {
            }
        }

        /// <summary>
        /// Processes a feature from a GTF gene model file.
        /// </summary>
        /// <param name="feature"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chrom"></param>
        /// <param name="attributes"></param>
        public void ProcessGtfFeature(MetadataListItem<List<string>> feature, long oneBasedStart, long oneBasedEnd, Chromosome chrom, Dictionary<string, string> attributes)
        {
            bool hasGeneId = attributes.TryGetValue("gene_id", out string geneId);
            bool hasTranscriptId = attributes.TryGetValue("transcript_id", out string transcriptId);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string proteinId);
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasSource = feature.SubItems.TryGetValue("source", out List<string> sourceish);
            bool hasStrand = feature.SubItems.TryGetValue("strand", out List<string> strandish);
            bool hasFrame = feature.SubItems.TryGetValue("frame", out List<string> framey);

            string source = hasSource ? sourceish[0] : "";
            if (!hasStrand) { return; } // strand is a required to do anything in this program
            string strand = strandish[0];
            int frame = 0;
            if (hasFrame) { int.TryParse(framey[0], out frame); }

            // Trim prefixes from the IDs
            string genePrefix = "gene:";
            string transcriptPrefix = "transcript:";
            if (hasGeneId && geneId.StartsWith(genePrefix))
            {
                string newGeneId = geneId.Substring(genePrefix.Length);
                feature.FreeText.Replace(geneId, newGeneId);
                geneId = newGeneId;
            }
            if (hasTranscriptId && transcriptId.StartsWith(transcriptPrefix))
            {
                string newTranscriptId = transcriptId.Substring(transcriptPrefix.Length);
                feature.FreeText.Replace(transcriptId, newTranscriptId);
                transcriptId = newTranscriptId;
            }
            if (hasProteinId && proteinId.StartsWith(transcriptPrefix))
            {
                proteinId = proteinId.Substring(transcriptPrefix.Length); // transcript id is used for protein id sometimes
            }

            // Catch the transcript features before they go by if available, i.e. if the file doesn't just have exons
            if (feature.Key == "transcript" && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chrom, source, strand, oneBasedStart, oneBasedEnd, feature);
                    Genes.Add(currentGene);
                    GenomeForest.Add(currentGene);
                }

                currentTranscript = new Transcript(transcriptId, currentGene, source, strand, oneBasedStart, oneBasedEnd, null, feature);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (feature.Key == "exon" || feature.Key == "CDS")
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chrom, source, strand, oneBasedStart, oneBasedEnd, feature);
                    Genes.Add(currentGene);
                    GenomeForest.Add(currentGene);
                }

                if (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID)
                {
                    if (currentTranscript != null)
                    {
                        Transcript.SetRegions(currentTranscript);
                        currentTranscript.FrameCorrection();
                    }
                    currentTranscript = new Transcript(transcriptId, currentGene, source, strand, oneBasedStart, oneBasedEnd, null, feature);
                    currentGene.Transcripts.Add(currentTranscript);
                    GenomeForest.Add(currentTranscript);
                }

                if (feature.Key == "exon")
                {
                    ISequence exon_dna = chrom.Sequence.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                    Exon exon = new Exon(currentTranscript, currentTranscript.IsStrandPlus() ? exon_dna : exon_dna.GetReverseComplementedSequence(),
                        source, oneBasedStart, oneBasedEnd, chrom.Sequence.ID, strand, feature);
                    if (exon.Length() > 0) { currentTranscript.Exons.Add(exon); }
                }
                else if (feature.Key == "CDS")
                {
                    CDS cds = new CDS(currentTranscript, chrom.Sequence.ID, source, strand, oneBasedStart, oneBasedEnd, frame);
                    if (hasProteinId) { currentTranscript.ProteinID = proteinId; }
                    if (cds.Length() > 0) { currentTranscript.CodingDomainSequences.Add(cds); }
                }
                else
                { // nothing to do
                }
            }
        }

        /// <summary>
        /// Splits attributes from free text, either from GTF or GFF
        /// </summary>
        /// <param name="freeText"></param>
        /// <returns></returns>
        public static Dictionary<string, string> SplitAttributes(string freeText)
        {
            Dictionary<string, string> attributes = new Dictionary<string, string>();
            foreach (string attrib in freeText.Split(';'))
            {
                string key;
                string val;
                if (freeText.Contains('=')) // GFF3
                {
                    key = attrib.Split('=')[0].TrimStart();
                    val = attrib.Split('=')[1].TrimStart();
                }
                else // GFF1 or GTF
                {
                    key = AttributeKey.Match(attrib.TrimStart()).Groups[1].Value;
                    val = AttributeValue.Match(attrib.TrimStart()).Groups[1].Value;
                }

                if (!attributes.TryGetValue(key, out string x)) // sometimes there are two tags, so avoid adding twice
                {
                    attributes.Add(key, val);
                }
            }
            return attributes;
        }

        /// <summary>
        /// Gets chromsomes as ISequence objects, containing feature metadata contained in the gene model
        /// </summary>
        /// <param name="geneModelFile"></param>
        /// <returns></returns>
        public static List<ISequence> SimplerParse(string geneModelFile)
        {
            ForceGffVersionTo2(geneModelFile, out string geneModelWithVersion2MarkedPath);
            List<ISequence> geneFeatures = new GffParser().Parse(geneModelWithVersion2MarkedPath).ToList();
            return geneFeatures;
        }

        /// <summary>
        /// Regex for getting the gff version
        /// </summary>
        private static Regex gffVersion = new Regex(@"(##gff-version\s+)(\d)");

        /// <summary>
        /// Required for using DotNetBio because it only handles GFF version 2 in the header.
        /// The only difference in the new version is within the attributes, which are stored as free text anyway.
        /// </summary>
        /// <param name="gffPath"></param>
        /// <param name="gffWithVersionMarked2Path"></param>
        public static void ForceGffVersionTo2(string gffPath, out string gffWithVersionMarked2Path)
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

        #endregion Methods -- Read Gene Model File

        #region Methods -- Applying Variants and Translation

        /// <summary>
        /// Create CDS for transcripts in this gene model, based on the translation from CDS in another model
        /// </summary>
        /// <param name="referenceGeneModel"></param>
        public void CreateCDSFromAnnotatedStartCodons(GeneModel referenceGeneModel)
        {
            referenceGeneModel.GenomeForest.Build(); // so we don't need to lock the IntervalTree if we end up parallelizing this method
            foreach (Gene g in Genes)
            {
                bool hasSource = g.FeatureMetadata.SubItems.TryGetValue("source", out List<string> sourceish);
                if (!referenceGeneModel.GenomeForest.Forest.TryGetValue(g.Chromosome.FriendlyName, out IntervalTree tree)) { continue; }
                foreach (Transcript t in g.Transcripts)
                {
                    List<Transcript> referenceTranscripts = tree.Query(t).OfType<Transcript>().ToList();
                    List<Transcript> referenceTranscriptsWithCDS = referenceTranscripts.Where(tt => tt.IsProteinCoding()).ToList();
                    foreach (Transcript tWithCds in referenceTranscriptsWithCDS)
                    {
                        lock (tWithCds)
                        {
                            if (t.CreateCDSFromAnnotatedStartCodons(tWithCds)) { break; } // for now, only use the first annotation found if any
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Print out to a
        /// </summary>
        /// <param name="outGffFilePath"></param>
        public void PrintToGTF(string outGffFilePath)
        {
            using (FileStream stream = new FileStream(outGffFilePath, FileMode.Create))
            {
                GffFormatter gff = new GffFormatter();
                foreach (Chromosome chrom in Genome.Chromosomes)
                {
                    if (GenomeForest.Forest.TryGetValue(chrom.FriendlyName, out var tree))
                    {
                        IEnumerable<Gene> genes = tree.Intervals.OfType<Gene>().OrderBy(g => g.OneBasedStart);
                        chrom.Sequence.Metadata["features"] = genes.SelectMany(g => g.GetFeatures()).ToList();
                        chrom.Sequence.ID = chrom.FriendlyName; // shortens to "1" from "1 dna:chromosome chromosome:GRCh37:1:1:249250621:1 REF"
                        gff.Format(stream, chrom.Sequence);
                    }
                }
            }
        }

        /// <summary>
        /// Creates UTRs for transcripts and intergenic regions after reading gene model
        /// </summary>
        public void CreateIntergenicRegions()
        {
            foreach (IntervalTree it in GenomeForest.Forest.Values)
            {
                Gene previousPositiveStrandGene = null;
                Gene previousNegativeStrandGene = null;

                // Create intergenic regions on each strand
                foreach (Gene gene in it.Intervals.OfType<Gene>().OrderBy(g => g.OneBasedStart))
                {
                    Intergenic intergenic = null;
                    Gene previous = gene.IsStrandPlus() ? previousPositiveStrandGene : previousNegativeStrandGene;
                    if (previous != null)
                    {
                        // if there's a previous gene, create the intergenic region
                        intergenic = new Intergenic(gene.Chromosome, gene.ChromosomeID, gene.Source, gene.Strand, previous.OneBasedEnd + 1, gene.OneBasedStart - 1);
                    }

                    // store previous genes on each strand
                    if (gene.IsStrandPlus())
                    {
                        previousPositiveStrandGene = gene;
                    }
                    if (gene.IsStrandMinus())
                    {
                        previousNegativeStrandGene = gene;
                    }

                    // add the intergenic region to the genome forest if it was created
                    if (intergenic != null && intergenic.Length() > 0)
                    {
                        GenomeForest.Add(intergenic);
                    }
                }
            }
        }

        public List<Protein> Translate(bool translateCodingDomains, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            return Genes.SelectMany(g => g.Translate(translateCodingDomains, incompleteTranscriptAccessions, selenocysteineContaining)).ToList();
        }

        /// <summary>
        /// Ensembl coding domain sequences (CDS) sometimes don't have start or stop codons annotated.
        /// The only way I can figure out how to tell which they are is to read in the protein FASTA and find the ones starting with X's or containing a stop codon '*'
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="proteinFastaPath"></param>
        /// <returns></returns>
        public static void GetImportantProteinAccessions(string proteinFastaPath, out Dictionary<string, string> proteinAccessionSequence, out HashSet<string> badProteinAccessions,
            out Dictionary<string, string> selenocysteineProteinAccessions)
        {
            Regex transcriptAccession = new Regex(@"(transcript:)([A-Za-z0-9_.]+)"); // need to include transcript accessions for when a GTF file is used and transcript IDs become the protein IDs
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(proteinFastaPath, true, DecoyType.None, false,
                ProteinDbLoader.EnsemblAccessionRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblFullNameRegex, ProteinDbLoader.EnsemblGeneNameRegex, null, out List<string> errors);
            proteinAccessionSequence = proteins.Select(p => new KeyValuePair<string, string>(p.Accession, p.BaseSequence))
                .Concat(proteins.Select(p => new KeyValuePair<string, string>(transcriptAccession.Match(p.FullName).Groups[2].Value, p.BaseSequence)))
                .ToDictionary(kv => kv.Key, kv => kv.Value);
            HashSet<string> badOnes = new HashSet<string>(proteins.Where(p => p.BaseSequence.Contains('X') || p.BaseSequence.Contains('*'))
                .SelectMany(p => new string[] { p.Accession, transcriptAccession.Match(p.FullName).Groups[2].Value }));
            badProteinAccessions = badOnes;
            selenocysteineProteinAccessions = proteins.Where(p => !badOnes.Contains(p.Accession) && p.BaseSequence.Contains('U')).ToDictionary(p => p.Accession, p => p.BaseSequence);
        }

        #endregion Methods -- Applying Variants and Translation
    }
}