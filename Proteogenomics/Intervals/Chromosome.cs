using Bio;
using System;

namespace Proteogenomics
{
    public class Chromosome
        : Interval
    {
        #region Public Properties

        public ISequence Sequence { get; set; }

        public string FriendlyName { get; set; }

        #endregion Public Properties

        #region Constructor

        public Chromosome(ISequence sequence)
            : base(sequence.ID, "+", 1, sequence.Count)
        {
            Sequence = sequence;
            FriendlyName = GetFriendlyChromosomeName(ChromosomeID);
        }

        #endregion Constructor

        #region Public Static Method

        public static string GetFriendlyChromosomeName(string chromosomeID)
        {
            return chromosomeID.Split(new string[] { ProteogenomicsUtility.ENSEMBL_FASTA_HEADER_DELIMETER }, StringSplitOptions.None)[0];
        }

        #endregion Public Static Method
    }
}