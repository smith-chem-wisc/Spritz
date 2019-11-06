using Bio;
using Bio.Algorithms.Translation;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class Translation
    {
        /// <summary>
        /// Stores accessions for checking that they are unique.
        /// </summary>
        private static HashSet<string> ProteinAccessions = new HashSet<string>();

        /// <summary>
        /// Translate the coding sequence of a transcript
        /// </summary>
        /// <param name="transcript"></param>
        /// <returns></returns>
        public static ISequence OneFrameTranslation(ISequence dnaSequence, bool mitochondrial)
        {
            ISequence rnaSequence = Transcription.Transcribe(dnaSequence);
            ISequence proteinSequence;
            if (!mitochondrial)
            {
                proteinSequence = ProteinTranslation.Translate(rnaSequence);
            }
            else
            {
                List<byte> aminoAcidSequence = new List<byte>();
                for (int codonNum = 0; codonNum < (int)(rnaSequence.Count / 3); codonNum++)
                {
                    // Check for start codon, and add 'M' if so
                    if (aminoAcidSequence.Count == 0 &&
                        CodonsVertebrateMitochondrial.START_CODONS.Contains(
                            new string(dnaSequence.GetSubSequence(codonNum * 3, GeneModel.CODON_SIZE).Select(bp => (char)bp).ToArray())))
                    {
                        aminoAcidSequence.Add((byte)CodonsVertebrateMitochondrial.DEFAULT_START_AA);
                        continue;
                    }
                    aminoAcidSequence.Add(CodonsVertebrateMitochondrial.TryLookup(rnaSequence, codonNum * 3, out byte aa) ? aa : (byte)'X');
                }
                proteinSequence = new Sequence(Alphabets.AmbiguousProtein, aminoAcidSequence.ToArray());
            }
            return proteinSequence;
        }

        public static string GetSafeProteinAccession(string accession)
        {
            string newAccession = accession;
            int i = 1;
            while (ProteinAccessions.Contains(newAccession))
            {
                newAccession = accession + "_" + i++.ToString();
            }
            ProteinAccessions.Add(newAccession);
            return newAccession;
        }

        /// <summary>
        /// Not used or tested right now...
        /// </summary>
        /// <param name="exons"></param>
        /// <param name="proteinID"></param>
        /// <returns></returns>
        public static Protein ThreeFrameTranslation(List<Exon> exons, string proteinID)
        {
            string seq = String.Join("", exons.Select(x => SequenceExtensions.ConvertToString(x.Sequence)));
            if (seq.Contains('N')) return null;
            ISequence dna_seq = new Sequence(Alphabets.DNA, seq);
            ISequence rna_seq = Transcription.Transcribe(exons[0].IsStrandPlus() ? dna_seq : dna_seq.GetReverseComplementedSequence());
            ISequence[] prot_seq = Enumerable.Range(0, 3).Select(i => ProteinTranslation.Translate(rna_seq, i)).ToArray();

            //return the protein sequence corresponding to the longest ORF
            return new Protein(prot_seq.SelectMany(s => SequenceExtensions.ConvertToString(s).Split('*')).OrderByDescending(s => s.Length).FirstOrDefault(), proteinID);
        }

        /// <summary>
        /// Get AA for a single DNA codon
        /// </summary>
        /// <param name="codon"></param>
        /// <param name="aa"></param>
        /// <returns></returns>
        public static bool TranslateDnaCodon(string codon, out byte aa)
        {
            ISequence rnaCodon = Transcription.Transcribe(new Sequence(Alphabets.DNA, codon.Select(c => (byte)c).ToArray()));
            return Codons.TryLookup(rnaCodon, 0, out aa);
        }
    }
}