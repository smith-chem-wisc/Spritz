namespace Genomics
{
    public class Chromosome
    {

        #region Public Constructors

        public string Name { get; set; }
        public NucleotideSequence Sequence { get; set; }
        public int Length { get { return Sequence.Length; } }

        #endregion Public Constructors

        #region Public Constructor

        public Chromosome(string name, char[] sequence)
        {
            this.Name = name;
            this.Sequence = new NucleotideSequence(sequence);
        } 

        #endregion

    }
}
