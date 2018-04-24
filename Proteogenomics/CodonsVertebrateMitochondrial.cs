using Bio;
using Bio.Algorithms.Translation;
using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public static class CodonsVertebrateMitochondrial
    {
        /// <summary>
        /// All start codons are translated as "M", even if alternative start codons, since an different tRNA is used. https://en.wikipedia.org/wiki/Start_codon
        /// </summary>
        public static readonly char DEFAULT_START_AA = 'M';

        /// <summary>
        /// Start codons for the vertebral mitochondrial genetic code.
        /// </summary>
        public static readonly HashSet<string> START_CODONS = new HashSet<string> { "ATT", "ATC", "ATA", "ATG", "GTG" };

        /// <summary>
        /// The mapping dictionary.
        /// </summary>
        private static readonly Dictionary<string, byte> CodonMap = new Dictionary<string, byte>();

        /// <summary>
        /// Lookup an amino acid based on a triplet of nucleotides. U U U for instance
        /// will result in Phenylalanine.  If the values cannot be
        /// found in the lookup table, <c>false</c> will be returned.
        /// </summary>
        /// <param name="n1">The first character.</param>
        /// <param name="n2">The second character.</param>
        /// <param name="n3">The third character.</param>
        /// <param name="aminoAcid">Mapped RNA value</param>
        /// <returns>True/False if the value exists</returns>
        public static bool TryLookup(byte n1, byte n2, byte n3, out byte aminoAcid)
        {
            var codon = new char[3];
            codon[0] = char.ToUpperInvariant((char)n1);
            codon[1] = char.ToUpperInvariant((char)n2);
            codon[2] = char.ToUpperInvariant((char)n3);
            return CodonMap.TryGetValue(new string(codon), out aminoAcid);
        }

        /// <summary>
        /// Lookup an amino acid within a sequence starting a certain offset.
        /// </summary>
        /// <param name="sequence">The source sequence to lookup from.</param>
        /// <param name="offset">
        /// The offset within the sequence from which to look at the next three
        /// nucleotides. Note that this offset begins its count at zero. Thus
        /// looking at a sequence described by "AUGGCG" and using a offset of 0
        /// would lookup the amino acid for codon "AUG" while using a offset of 1
        /// would lookup the amino acid for codon "UGG".
        /// </param>
        /// <returns>An amino acid from the protein alphabet</returns>
        public static byte Lookup(ISequence sequence, int offset)
        {
            byte value;
            if (TryLookup(sequence, offset, out value))
            {
                return value;
            }

            throw new InvalidOperationException(
                string.Format(null, "({0},{1},{2}) not found.", (char)sequence[offset], (char)sequence[offset + 1], (char)sequence[offset + 2]));
        }

        /// <summary>
        /// Tries to lookup an amino acid within a sequence starting a certain offset.
        /// </summary>
        /// <param name="sequence">The source sequence to lookup from.</param>
        /// <param name="offset">
        /// The offset within the sequence from which to look at the next three
        /// nucleotides. Note that this offset begins its count at zero. Thus
        /// looking at a sequence described by "AUGGCG" and using a offset of 0
        /// would lookup the amino acid for codon "AUG" while using a offset of 1
        /// would lookup the amino acid for codon "UGG".
        /// </param>
        /// <param name="aminoAcid">An amino acid from the protein alphabet</param>
        /// <returns><c>true</c>, if the triplet of nucleotides could
        /// be mapped, <c>false</c> otherwise</returns>
        public static bool TryLookup(ISequence sequence, int offset, out byte aminoAcid)
        {
            if (sequence == null)
            {
                throw new ArgumentNullException("sequence");
            }
            if (offset >= sequence.Count - 2)
            {
                throw new ArgumentException("offset overflow");
            }
            return TryLookup(sequence[offset], sequence[offset + 1], sequence[offset + 2], out aminoAcid);
        }

        #region Constructor

        /// <summary>
        /// Initializes the Codon map dictionary.
        /// </summary>
        static CodonsVertebrateMitochondrial()
        {
            string[] codons = "TTT/F,TTC/F,TTA/L,TTG/L,TCT/S,TCC/S,TCA/S,TCG/S,TAT/Y,TAC/Y,TAA/*,TAG/*,TGT/C,TGC/C,TGA/W,TGG/W,CTT/L,CTC/L,CTA/L,CTG/L,CCT/P,CCC/P,CCA/P,CCG/P,CAT/H,CAC/H,CAA/Q,CAG/Q,CGT/R,CGC/R,CGA/R,CGG/R,ATT/I+,ATC/I+,ATA/M+,ATG/M+,ACT/T,ACC/T,ACA/T,ACG/T,AAT/N,AAC/N,AAA/K,AAG/K,AGT/S,AGC/S,AGA/*,AGG/*,GTT/V,GTC/V,GTA/V,GTG/V+,GCT/A,GCC/A,GCA/A,GCG/A,GAT/D,GAC/D,GAA/E,GAG/E,GGT/G,GGC/G,GGA/G,GGG/G"
                .Split(',');

            foreach (string codon in codons)
            {
                string[] codonAndAA = codon.Split('/');
                string rnaCodon = Transcription.Transcribe(new Sequence(Alphabets.DNA, codonAndAA[0])).ToString();
                CodonMap.Add(rnaCodon, (byte)codonAndAA[1][0]);
            }

            CodonMap.Add("---", Alphabets.Protein.Gap);
        }

        #endregion Constructor
    }
}