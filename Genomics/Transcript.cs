using Bio;
using Bio.Algorithms.Translation;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class Transcript
    {

        #region Private Properties

        private MetadataListItem<List<string>> metadata { get; set; }

        #endregion

        #region Public Properties

        public string ID { get; set; }
        public string Strand { get; set; }
        public Gene Gene { get; set; }
        public List<Exon> Exons { get; set; } = new List<Exon>();
        public List<Exon> CDS { get; set; } = new List<Exon>();
        public string ProteinID { get; set; }

        #endregion Public Properties

        #region Public Constructors

        public Transcript(string ID, Gene gene, MetadataListItem<List<string>> metadata, string ProteinID = null)
        {
            this.ID = ID;
            this.ProteinID = ProteinID == null ? ID : ProteinID;
            this.Gene = gene;
            this.metadata = metadata;
            this.Strand = metadata.SubItems["strand"][0];
        }

        #endregion Public Constructors

        #region Public Methods

        public Protein translate(bool translateCDS, bool includeVariants)
        {
            if (ProteinID == "ENSP00000473460")
                ProteinID = ProteinID;
            return one_frame_translation(translateCDS, String.Join("", (translateCDS ? CDS : Exons).Select(x => SequenceExtensions.ConvertToString(x.Sequence))));
        }

        private static Dictionary<Tuple<string, string>, List<Exon>> chromID_strand_binnedCDS = new Dictionary<Tuple<string, string>, List<Exon>>();
        public static Dictionary<Tuple<string, string>, List<Exon>> binCDS_by_chromID_strand(GeneModel geneModelWithCDS)
        {
            Dictionary<Tuple<string, string>, List<Exon>> binnedCDS = new Dictionary<Tuple<string, string>, List<Exon>>();
            foreach (Exon x in geneModelWithCDS.genes.SelectMany(g => g.transcripts).SelectMany(t => t.CDS))
            {
                Tuple<string, string> key = new Tuple<string, string>(x.ChromID, x.Strand);
                binnedCDS.TryGetValue(key, out List<Exon> xs);
                if (xs == null) binnedCDS.Add(key, new List<Exon> { x });
                else binnedCDS[key].Add(x);
            }
            return binnedCDS;
        }

        public Protein translate(Dictionary<Tuple<string, string>, List<Exon>> chromID_strand_binnedCDS, int min_length, bool includeVariants)
        {
            List<Exon> annotated_starts = chromID_strand_binnedCDS[new Tuple<string,string>(Gene.ChromID, Strand)].Where(x => 
                Exons.Any(xx => 
                    xx.includes(Strand == "+" ? x.OneBasedStart : x.OneBasedEnd) // must include the start of the stop codon
                    && xx.includes(Strand == "+" ? x.OneBasedStart + 2 : x.OneBasedEnd - 2))).ToList(); // and the end of the stop codon
            char terminating_character = ProteinAlphabet.Instance.GetFriendlyName(Alphabets.Protein.Ter)[0];
            if (annotated_starts.Count > 0)
            {
                // gets the first annotated start that produces
                foreach (Exon annotated_start in annotated_starts)
                {
                    long start_codon_start = Strand == "+" ? annotated_start.OneBasedStart : annotated_start.OneBasedEnd; // CDS on the reverse strand have start and end switched
                    string seq = String.Join("", Exons.Select(x => SequenceExtensions.ConvertToString(x.Sequence)));
                    if (Strand == "+")
                    {
                        long exon_length_before = Exons.Where(x => x.OneBasedEnd < annotated_start.OneBasedStart).Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                        long exon_zero_based_start = start_codon_start - Exons.FirstOrDefault(x => x.includes(annotated_start.OneBasedStart)).OneBasedStart;
                        long ZeroBasedStart = exon_length_before + exon_zero_based_start;
                        long length_after = seq.Length - ZeroBasedStart;
                        if (seq.Contains('N')) return null;
                        string subseq = SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, seq).GetSubSequence(ZeroBasedStart, length_after));
                        Protein p = one_frame_translation(false, subseq);
                        if (p.BaseSequence.Length >= min_length) return p;
                        else continue;
                    }
                    else
                    {
                        long length = Exons.Sum(x => x.OneBasedEnd - x.OneBasedStart + 1);
                        long chop = Exons.Where(x => x.OneBasedEnd >= annotated_start.OneBasedEnd).Sum(x => annotated_start.OneBasedEnd < x.OneBasedStart ? x.OneBasedEnd - x.OneBasedStart + 1 : x.OneBasedEnd - annotated_start.OneBasedEnd);
                        if (seq.Contains('N')) return null;
                        string subseq = SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, seq).GetSubSequence(0, length - chop));
                        Protein p = one_frame_translation(false, subseq);
                        if (p.BaseSequence.Length >= min_length) return p;
                        else continue;
                    }
                }
            }
            return null;
            //return three_frame_translation();
        }

        #endregion Public Methods

        #region Private Methods

        private Protein one_frame_translation(bool translateCDS, string seq)
        {
            if (seq.Length == 0 || !translateCDS && Exons.Count == 0 || translateCDS && CDS.Count == 0 || seq.Contains('N'))
                return null;
            ISequence dna_seq = new Sequence(Alphabets.DNA, seq);
            ISequence rna_seq = Transcription.Transcribe((translateCDS ? CDS : Exons)[0].Strand == "+" ? dna_seq : dna_seq.GetReverseComplementedSequence());
            ISequence prot_seq = ProteinTranslation.Translate(rna_seq);
            return new Protein(SequenceExtensions.ConvertToString(prot_seq).Split('*')[0], ProteinID);
        }

        private Protein three_frame_translation()
        {
            string seq = String.Join("", Exons.Select(x => SequenceExtensions.ConvertToString(x.Sequence)));
            if (seq.Contains('N')) return null;
            ISequence dna_seq = new Sequence(Alphabets.DNA, seq);
            ISequence rna_seq = Transcription.Transcribe(Exons[0].Strand == "+" ? dna_seq : dna_seq.GetReverseComplementedSequence());
            ISequence[] prot_seq = Enumerable.Range(0, 3).Select(i => ProteinTranslation.Translate(rna_seq, i)).ToArray();

            //return the protein sequence corresponding to the longest ORF
            return new Protein(prot_seq.SelectMany(s => SequenceExtensions.ConvertToString(s).Split('*')).OrderByDescending(s => s.Length).FirstOrDefault(), ProteinID);
        }

        #endregion Private Methods

    }
}
