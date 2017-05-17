using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class SequenceVariant : ChromosomeSegment
    {

        #region Public Properties

        /// <summary>
        /// Phred score for this variant
        /// </summary>
        public double Qual { get; set; }

        /// <summary>
        /// Reference Allele
        /// </summary>
        public string Ref { get; set; }

        /// <summary>
        /// Alternate allele or comma-separated alleles
        /// </summary>
        public string Alt { get; set; }

        public string Filter { get; set; }

        public Dictionary<string, string> info { get; set; }

        public new int length { get { return Math.Max(Ref.Length, Alt.Length); } }

        #endregion Public Properties

        #region Public Constructor

        public SequenceVariant(Chromosome chrom, int position, string id, string reference, string alternate, double qual, string filter, Dictionary<string, string> info)
            : base(id, chrom, "+", position, position + Math.Max(reference.Length, alternate.Length) - 1, null, null)
        {
            this.Qual = qual;
            this.Ref = reference;
            this.Alt = alternate;
            this.info = info;
        }

        #endregion Public Constructor

        #region Public Methods

        public override string ToString()
        {
            return String.Join("\t", new string[] { chrom.name, start.ToString(), id, Ref, Alt, Qual.ToString(), Filter.ToString(), String.Join(";", info.Select(x => x.Key + "=" + x.Value)) });
        }

        #endregion Public Methods

    }
}
