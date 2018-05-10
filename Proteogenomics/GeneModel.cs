using Bio;
using Bio.IO.Gff;
using Bio.VCF;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Bio.Extensions;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace Proteogenomics
{
    /// <summary>
    /// Contains representation of genes, transcripts, exons, etc. represented in a gene model. Can be amended with variants.
    /// </summary>
    public class GeneModel
    {
        /// <summary>
        /// Gets the first instance of a word
        /// </summary>
        private static Regex AttributeKey = new Regex(@"([\w]+)");

        /// <summary>
        /// Gets anything inside quotes
        /// </summary>
        private static Regex AttributeValue = new Regex(@"""([\w.]+)""");

        /// <summary>
        /// Used to check if variants were added before applying them
        /// </summary>
        private static List<Variant> PreviouslyAddedVariants;

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

        private Gene currentGene = null;
        private Transcript currentTranscript = null;

        public void ReadGeneFeatures(string geneModelFile)
        {
            ForceGffVersionTo2(geneModelFile, out string geneModelWithVersion2MarkedPath);
            List<ISequence> geneFeatures = new GffParser().Parse(geneModelWithVersion2MarkedPath).ToList();

            foreach (ISequence chromFeatures in geneFeatures)
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
                            key = AttributeKey.Match(attrib.TrimStart()).Groups[1].Value;
                            val = AttributeValue.Match(attrib.TrimStart()).Groups[1].Value;
                        }

                        if (!attributes.TryGetValue(key, out string x)) // sometimes there are two tags, so avoid adding twice
                        {
                            attributes.Add(key, val);
                        }
                    }

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
            bool hasTranscriptVersion = attributes.TryGetValue("version", out string transcriptVersion) && hasTranscriptId;
            bool hasExonId = attributes.TryGetValue("exon_id", out string exonId);
            bool hasProteinId = attributes.TryGetValue("protein_id", out string proteinId);
            bool hasStrand = feature.SubItems.TryGetValue("strand", out List<string> strandish); // false if empty ("." in GFF format)
            bool hasFrame = feature.SubItems.TryGetValue("frame", out List<string> framey); // false if empty ("." in GFF format)
            if (!hasStrand) { return; } // strand is a required to do anything in this program
            string strand = strandish[0];
            int frame = 0;
            if (hasFrame) { int.TryParse(framey[0], out frame); }

            if (hasGeneId && (currentGene == null || hasGeneId && geneId != currentGene.ID))
            {
                currentGene = new Gene(geneId, chrom, strand, oneBasedStart, oneBasedEnd);
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
                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd, null, null);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (hasExonId)
            {
                ISequence exon_dna = chrom.Sequence.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                Exon exon = new Exon(currentTranscript, currentTranscript.IsStrandPlus() ? exon_dna : exon_dna.GetReverseComplementedSequence(),
                    oneBasedStart, oneBasedEnd, chrom == null ? "" : chrom.ChromosomeID, strand, null);
                if (exon.Length() > 0)
                {
                    currentTranscript.Exons.Add(exon);
                }
            }
            else if (hasProteinId)
            {
                CDS cds = new CDS(currentTranscript, chrom.Sequence.ID, strand, oneBasedStart, oneBasedEnd, null, frame);
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
            bool hasStrand = feature.SubItems.TryGetValue("strand", out List<string> strandish);
            bool hasFrame = feature.SubItems.TryGetValue("frame", out List<string> framey);
            if (!hasStrand) { return; } // strand is a required to do anything in this program
            string strand = strandish[0];
            int frame = 0;
            if (hasFrame) { int.TryParse(framey[0], out frame); }

            // Catch the transcript features before they go by if available, i.e. if the file doesn't just have exons
            if (feature.Key == "transcript" && (currentTranscript == null || hasTranscriptId && transcriptId != currentTranscript.ID))
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chrom, strand, oneBasedStart, oneBasedEnd);
                    Genes.Add(currentGene);
                    GenomeForest.Add(currentGene);
                }

                currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd, null, null);
                currentGene.Transcripts.Add(currentTranscript);
                GenomeForest.Add(currentTranscript);
            }

            if (feature.Key == "exon" || feature.Key == "CDS")
            {
                if (currentGene == null || hasGeneId && geneId != currentGene.ID)
                {
                    currentGene = new Gene(geneId, chrom, strand, oneBasedStart, oneBasedEnd);
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
                    currentTranscript = new Transcript(transcriptId, transcriptVersion, currentGene, strand, oneBasedStart, oneBasedEnd, null, null);
                    currentGene.Transcripts.Add(currentTranscript);
                    GenomeForest.Add(currentTranscript);
                }

                if (feature.Key == "exon")
                {
                    ISequence exon_dna = chrom.Sequence.GetSubSequence(oneBasedStart - 1, oneBasedEnd - oneBasedStart + 1);
                    Exon exon = new Exon(currentTranscript, currentTranscript.IsStrandPlus() ? exon_dna : exon_dna.GetReverseComplementedSequence(),
                        oneBasedStart, oneBasedEnd, chrom.Sequence.ID, strand, null);
                    if (exon.Length() > 0)
                    {
                        currentTranscript.Exons.Add(exon);
                    }
                }
                else if (feature.Key == "CDS")
                {
                    CDS cds = new CDS(currentTranscript, chrom.Sequence.ID, strand, oneBasedStart, oneBasedEnd, null, frame);
                    if (cds.Length() > 0)
                    {
                        currentTranscript.CodingDomainSequences.Add(cds);
                    }
                }
                else
                { // nothing to do
                }
            }
        }

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

        #endregion Methods -- Read Gene Model File

        #region Methods -- Applying Variants and Translation

        /// <summary>
        /// Create CDS for transcripts in this gene model, based on the translation from CDS in another model
        /// </summary>
        /// <param name="referenceGeneModel"></param>
        public void CreateCDSFromAnnotatedStartCodons(GeneModel referenceGeneModel)
        {
            foreach (Gene g in Genes)
            {
                if (!referenceGeneModel.GenomeForest.Forest.TryGetValue(g.Chromosome.FriendlyName, out IntervalTree tree)) { continue; }
                foreach (Transcript t in g.Transcripts)
                {
                    List<Transcript> referenceTranscripts = tree.Query(t).OfType<Transcript>().ToList();
                    List<Transcript> referenceTranscriptsWithCDS = referenceTranscripts.Where(tt => tt.IsProteinCoding()).ToList();
                    t.CreateCDSFromAnnotatedStartCodons(referenceTranscriptsWithCDS.FirstOrDefault());
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
                        intergenic = new Intergenic(gene.Chromosome, gene.ChromosomeID, gene.Strand, previous.OneBasedEnd + 1, gene.OneBasedStart - 1, null);
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

        /// <summary>
        /// Add variants to relevant genomic regions and annotate them if on transcripts
        /// </summary>
        /// <param name="variants"></param>
        public void AddVariantAnnotations(List<Variant> variants)
        {
            Parallel.ForEach(variants.OrderByDescending(v => v.OneBasedStart).ToList(), v =>
            {
                if (GenomeForest.Forest.TryGetValue(Chromosome.GetFriendlyChromosomeName(v.ChromosomeID), out var intervalTree))
                {
                    foreach (Interval i in intervalTree.Stab(v.OneBasedStart))
                    {
                        lock (i)
                        {
                            i.Variants.Add(v);
                            Transcript t = i as Transcript;
                            if (t != null)
                            {
                                t.AnnotateWithVariant(v);
                            }
                        }
                    }
                }
            });
            PreviouslyAddedVariants = variants;
        }

        /// <summary>
        /// Apply a list of variants to intervals within this gene model
        /// </summary>
        /// <param name="variants"></param>
        /// <returns></returns>
        public List<Transcript> ApplyVariants(List<Variant> variants)
        {
            // Must add variants before applying them to this gene model
            if (PreviouslyAddedVariants != variants) {  AddVariantAnnotations(variants); }
            return ApplyVariants(variants, Genes.SelectMany(g => g.Transcripts).ToList());
        }

        /// <summary>
        /// Apply variants to transcripts
        /// </summary>
        /// <param name="variants"></param>
        public static List<Transcript> ApplyVariants(List<Variant> variants, List<Transcript> transcripts)
        {
            return transcripts.SelectMany(t => ApplyVariantsCombinitorially(t)).ToList();
        }

        /// <summary>
        /// Apply variants to a transcript
        /// </summary>
        /// <param name="t"></param>
        /// <returns></returns>
        public static List<Transcript> ApplyVariantsCombinitorially(Transcript t)
        {
            // Clear out annotations from non-combinitoric add/annotate method
            t.VariantAnnotations.Clear();
            t.ProteinSequenceVariations.Clear();

            List<Transcript> newTranscripts = new List<Transcript> { new Transcript(t) };
            int heterozygousNonsynonymousCount = 0;
            List<Variant> transcriptVariants = t.Variants.OrderByDescending(v => v.OneBasedStart).ToList(); // reversed, so that the coordinates of each successive variant is not changed
            foreach (Variant v in transcriptVariants)
            {
                List<Transcript> newOnes = new List<Transcript>();
                foreach (var nt in newTranscripts)
                {
                    var newerOnes = nt.ApplyVariantCombinitorics(v, out var effects); // expands only when there is a heterozygous nonsynonymous variation
                    bool heterozygousNonsynonymous = v.GenotypeType == GenotypeType.HETEROZYGOUS && effects.Effects.Any(eff => eff.IsNonsynonymous());
                    if (heterozygousNonsynonymous)
                    {
                        heterozygousNonsynonymousCount++;
                    }
                    if (heterozygousNonsynonymousCount > 5 && heterozygousNonsynonymous)
                    {
                        Transcript.combinatoricFailures.Add("Heterozygous nonsynonymous variant: transcript:" + v.ToString() + " wasn't included in transcript: " + t.ID + " or protein: " + t.ProteinID);
                        break; // avoid large combinitoric problems for now (heterozygous, nonsynonymous count > 5), but still stick in the homozygous variation
                    }
                    newOnes.AddRange(newerOnes);
                }
                newTranscripts = newOnes;
            }
            return newTranscripts;
        }

        public List<Protein> Translate(bool translateCodingDomains, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            return Genes.SelectMany(g => g.Translate(translateCodingDomains, incompleteTranscriptAccessions, selenocysteineContaining)).ToList();
        }

        /// <summary>
        /// Merges another gene model into this one
        /// </summary>
        /// <param name="model"></param>
        public void Merge(GeneModel model)
        {
            // Get the nearby intervals based on the interval forest before changing the forest
            List<Tuple<Gene, List<Interval>>> nearbyReferenceIntervals = new List<Tuple<Gene, List<Interval>>>();
            foreach (Gene g in model.Genes)
            {
                if (GenomeForest.Forest.TryGetValue(g.Chromosome.FriendlyName, out var tree))
                {
                    var nearbyIntervals = tree.Query(g).ToList();
                    nearbyReferenceIntervals.Add(new Tuple<Gene, List<Interval>>(g, nearbyIntervals));
                }
            }

            // Iterate over the possible new genes
            foreach (var geneAndNearby in nearbyReferenceIntervals)
            {
                Gene newGene = geneAndNearby.Item1;
                List<Gene> nearbyGenes = geneAndNearby.Item2.OfType<Gene>().ToList();
                List<Transcript> nearbyTranscripts = geneAndNearby.Item2.OfType<Transcript>().ToList();

                // If there are no nearby genes, add the gene
                if (nearbyGenes.Count == 0)
                {
                    Genes.Add(newGene);
                    GenomeForest.Add(newGene);
                    foreach (var t in newGene.Transcripts)
                    {
                        GenomeForest.Add(t);
                    }
                    continue;
                }

                // Otherwise, check if the transcript is already in this reference (same # exons and mRNA sequence)
                // Then, find the most similar gene, i.e. having the least overlap with this possible new gene
                foreach (Transcript t in newGene.Transcripts)
                {
                    bool sameAsReference = nearbyTranscripts.Any(refT => refT.Exons.Count == t.Exons.Count && refT.SplicedRNA().ConvertToString() == t.SplicedRNA().ConvertToString());
                    if (sameAsReference) { continue; }
                    Gene mostSimilar = nearbyGenes.OrderBy(g => g.Minus(newGene).Sum(outside => outside.Length())).First();
                    mostSimilar.Transcripts.Add(t);
                    GenomeForest.Add(t);
                }
            }
        }

        #endregion Methods -- Applying Variants and Translation
    }
}