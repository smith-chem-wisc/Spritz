using Bio.Algorithms.Translation;
using System;

namespace Proteogenomics
{
    public static class CodonExtensions
    {
        /// <summary>
        /// Translate a codon into an amino acid
        /// </summary>
        /// <param name="mitochondrial"></param>
        /// <param name="codon"></param>
        /// <param name="aminoAcid"></param>
        /// <returns></returns>
        public static bool TryTranslateCodon(bool mitochondrial, string codon, out byte aminoAcid)
        {
            if (codon.Length != GeneModel.CODON_SIZE)
            {
                throw new ArgumentException("Codon size not supported: " + codon);
            }
            return TryTranslateBytes(mitochondrial, (byte)codon[0], (byte)codon[1], (byte)codon[2], out aminoAcid);
        }

        /// <summary>
        /// Translate 3 byte representations of DNA bases to an amino acid
        /// </summary>
        /// <param name="mitochondrial"></param>
        /// <param name="base1"></param>ec
        /// <param name="base2"></param>
        /// <param name="base3"></param>
        /// <param name="aminoAcid"></param>
        /// <returns></returns>
        public static bool TryTranslateBytes(bool mitochondrial, byte base1, byte base2, byte base3, out byte aminoAcid)
        {
            if (mitochondrial)
            {
                return CodonsVertebrateMitochondrial.TryLookup(
                    Transcription.GetRnaComplement(base1),
                    Transcription.GetRnaComplement(base2),
                    Transcription.GetRnaComplement(base3),
                    out aminoAcid);
            }
            else
            {
                return Codons.TryLookup(
                    Transcription.GetRnaComplement(base1),
                    Transcription.GetRnaComplement(base2),
                    Transcription.GetRnaComplement(base3),
                    out aminoAcid);
            }
        }
    }
}