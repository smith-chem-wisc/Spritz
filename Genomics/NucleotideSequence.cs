using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class NucleotideSequence
    {

        #region Public Properties

        public char[] Sequence { get; set; }
        //public int OneBasedStart { get; set; } = -1;
        //public int OneBasedEnd { get; set; } = -1;
        //public int start_codon_onebasedstart { get; set; } = -1;
        //public int stop_codon_onebasedstart { get; set; } = -1;
        public int Length { get { return Sequence.Length; } }

        #endregion Public Properties

        #region Public Constructor

        public NucleotideSequence(char[] sequence)
        {
            this.Sequence = sequence;
        }

        public NucleotideSequence(string sequence)
        {
            this.Sequence = sequence.ToCharArray();
        }

        #endregion Public Constructor

        #region IUPAC Codes

        public static char[] amino_acids = "ACDEFGHIKLMNPQRSTVWY".ToCharArray();

        public static Dictionary<char, string> amino_acids_1to3 = new Dictionary<char, string>()
        {
            { 'A', "Ala" },
            { 'C', "Cys" },
            { 'D', "Asp" },
            { 'E', "Glu" },
            { 'F', "Phe" },
            { 'G', "Gly" },
            { 'H', "His" },
            { 'I', "Ile" },
            { 'K', "Lys" },
            { 'L', "Leu" },
            { 'M', "Met" },
            { 'N', "Asn" },
            { 'P', "Pro" },
            { 'Q', "Gln" },
            { 'R', "Arg" },
            { 'S', "Ser" },
            { 'T', "Thr" },
            { 'V', "Val" },
            { 'W', "Trp" },
            { 'Y', "Tyr" }
        };
        public static Dictionary<string, char> amino_acids_3to1 =
            amino_acids_1to3.ToList().ToDictionary(
                one_to_three => one_to_three.Value,
                one_to_three => one_to_three.Key);


        public static char[] ambiguous_dna_letters = "GATCRYWSMKHBVDN".ToCharArray();
        public static char[] unambiguous_dna_letters = "GATC".ToCharArray();
        public static Dictionary<char, string> ambiguous_dna_values = new Dictionary<char, string>()
        {
            { 'A', "A" },
            { 'C', "C" },
            { 'G', "G" },
            { 'T', "T" },
            { 'M', "AC" },
            { 'R', "AG" },
            { 'W', "AT" },
            { 'S', "CG" },
            { 'Y', "CT" },
            { 'K', "GT" },
            { 'V', "ACG" },
            { 'H', "ACT" },
            { 'D', "AGT" },
            { 'B', "CGT" },
            { 'X', "GATC" },
            { 'N', "GATC" }
        };

        public static Dictionary<char, char> ambiguous_dna_complement = new Dictionary<char, char>()
        {
            { 'A', 'T' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'T', 'A' },
            { 'M', 'K' },
            { 'R', 'Y' },
            { 'W', 'W' },
            { 'S', 'S' },
            { 'Y', 'R' },
            { 'K', 'M' },
            { 'V', 'B' },
            { 'H', 'D' },
            { 'D', 'H' },
            { 'B', 'V' },
            { 'X', 'X' },
            { 'N', 'N' }
        };

        public static Dictionary<char, char> unambiguous_dna_complement = new Dictionary<char, char>()
        {
            { 'A', 'T' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'T', 'A' },
        };

        public static List<string> standard_start_codons = new List<string>() { "ATG" };
        public static List<string> extended_start_codons = new List<string>() { "TTG", "CTG", "ATG" };
        public static List<string> standard_stop_codons = new List<string>() { "TAA", "TAG", "TGA" };
        public static List<string> mitochon_start_codons = new List<string>() { "ATT", "ATC", "ATA", "ATG", "GTG" };
        public static List<string> mitochon_stop_codons = new List<string>() { "TAA", "TAG", "AGA", "AGG" };

        public static char[] standard_base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG".ToCharArray();
        public static char[] standard_base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG".ToCharArray();
        public static char[] standard_base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG".ToCharArray();
        public static char[] standard_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".ToCharArray();
        public static char[] mitochon_acids = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG".ToCharArray();

        public static Dictionary<string, char> standard_code =
            Enumerable.Range(0, standard_acids.Length).ToDictionary(
                i => standard_base1[i].ToString() + standard_base2[i].ToString() + standard_base3[i].ToString(),
                i => standard_acids[i]);

        public static Dictionary<string, char> mitochon_code =
            Enumerable.Range(0, mitochon_acids.Length).ToDictionary(
                i => standard_base1[i].ToString() + standard_base2[i].ToString() + standard_base3[i].ToString(),
                i => mitochon_acids[i]);

        #endregion

        #region Public Methods

        public string get_range(int oneBasedStart, int oneBasedEnd)
        {
            if (oneBasedStart > oneBasedEnd) // this is on the reverse strand
            {
                int start = oneBasedStart;
                oneBasedStart = oneBasedEnd;
                oneBasedEnd = start;
            }
            else if (oneBasedStart == oneBasedEnd)
            {
                return "";
            }
            char[] result_arr = new char[oneBasedEnd - oneBasedStart + 1];
            Array.Copy(Sequence, oneBasedStart - 1, result_arr, 0, oneBasedEnd - oneBasedStart + 1);
            return new string(result_arr);
        }

        public Protein get_peptide(Dictionary<string, char> genetic_code, string id)
        {
            char[] result_arr = new char[Length / 3];
            for (int codon_idx = 0; codon_idx < Length / 3; codon_idx++)
            {
                char[] codon = new char[3];
                Array.Copy(Sequence, codon_idx * 3, codon, 0, 3);
                result_arr[codon_idx] = genetic_code[new string(codon)];
            }
            return new Protein(new string(result_arr), id);
        }

        #endregion Public Methods

    }
}
