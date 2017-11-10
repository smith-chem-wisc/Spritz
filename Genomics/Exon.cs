using Bio;
using Bio.VCF;
using Bio.Extensions;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class Exon
    {

        #region Private Properties

        private MetadataListItem<List<string>> Metadata { get; set; }

        #endregion Private Properties

        #region Public Properties

        public ISequence Sequence { get; set; }

        public string ChromID { get; set; }

        public string Strand { get; set; }

        public long OneBasedStart { get; set; } = -1;

        public long OneBasedEnd { get; set; } = -1;

        /// <summary>
        /// Used for both objects straight from VCF and for specific Variants, which extend VariantContext
        /// </summary>
        public List<VariantContext> Variants { get; set; } = new List<VariantContext>();

        #endregion Public Properties

        #region Public Constructor

        public Exon(ISequence Sequence, long start, long stop, string chrom_id, MetadataListItem<List<string>> metadata)
        {
            this.Sequence = Sequence;
            this.OneBasedStart = start;
            this.OneBasedEnd = stop;
            this.ChromID = chrom_id;
            this.Strand = metadata.SubItems["strand"][0];
            this.Metadata = metadata;
        }

        #endregion Public Constructor

        #region Private Constructor

        private Exon(Exon exon)
        {
            this.Sequence = exon.Sequence;
            this.OneBasedStart = exon.OneBasedStart;
            this.OneBasedEnd = exon.OneBasedEnd;
            this.ChromID = exon.ChromID;
            this.Strand = exon.Metadata.SubItems["strand"][0];
            this.Metadata = exon.Metadata;
        }

        #endregion Private Constructor

        #region Public Methods

        public IEnumerable<Exon> GetExonSequences(bool getVariantSequences = false, double homozygousThreshold = 1)
        {
            string seq = SequenceExtensions.ConvertToString(this.Sequence);
            if (Variants.Count == 0 || !getVariantSequences)
                return new List<Exon> { new Exon(this) };

            List<Variant> homozygous = new List<Variant>();
            List<Variant> heterozygous = new List<Variant>();
            foreach (Variant var in Variants.SelectMany(v => Variant.ParseVariantContext(v)))
            {
                // todo: allow depth filter here using "DP" column, especially for indels
                if (var.AlleleFrequency >= homozygousThreshold) homozygous.Add(var); // homozygous if the allele frequency is about 1
                else heterozygous.Add(var);
            }

            List<List<Variant>> possibleHaplotypes = new List<List<Variant>> { new List<Variant>(homozygous) };
            possibleHaplotypes.AddRange(
                Enumerable.Range(1, heterozygous.Count).SelectMany(k =>
                    ExtensionMethods.Combinations(heterozygous, k)
                        .Select(heteroAlleles => new List<Variant>(homozygous).Concat(heteroAlleles).ToList()))
                    .ToList());

            List<Exon> sequences = new List<Exon>();
            foreach (List<Variant> haplotype in possibleHaplotypes)
            {
                long tmpOneBasedStart = OneBasedStart;
                StringBuilder builder = new StringBuilder();
                haplotype.OrderBy(v => v.OneBasedStart);

                int subseqStart;
                int subseqCount;
                foreach (Variant var in haplotype)
                {
                    if (var.OneBasedStart < tmpOneBasedStart)
                        continue; // variant collides with previous variant

                    subseqStart = (int)(tmpOneBasedStart - this.OneBasedStart);
                    subseqCount = (int)(var.OneBasedStart - tmpOneBasedStart);
                    if (subseqStart + subseqCount - 1 > seq.Length) // allele extends beyond exon, so mind the exon boundary
                    {
                        builder.Append(seq.Substring(subseqStart));
                        builder.Append(var.AlternateAllele.Substring(0, (int)(this.OneBasedEnd - var.OneBasedStart)));
                    }
                    else
                    {
                        builder.Append(seq.Substring(subseqStart, subseqCount));
                        builder.Append(var.AlternateAllele);
                    }
                    tmpOneBasedStart = var.OneBasedStart + var.ReferenceAllele.Length;
                }

                // Append the remaining sequence
                subseqStart = (int)(tmpOneBasedStart - this.OneBasedStart);
                subseqCount = (int)(OneBasedEnd - tmpOneBasedStart + 1);
                builder.Append(seq.Substring(subseqStart, subseqCount));
                Exon x = new Exon(this);
                x.Sequence = new Sequence(Alphabets.DNA, builder.ToString());
                x.Variants = haplotype.OfType<VariantContext>().ToList();
                sequences.Add(x);
            }
            return sequences;
        }

        //public override int GetHashCode()
        //{
        //    return metadata.Key ^ metadata.SubItems["features"];
        //}

        public bool IsBefore(Exon segment)
        {
            return OneBasedStart < segment.OneBasedStart && OneBasedEnd < segment.OneBasedEnd && OneBasedEnd < segment.OneBasedStart;
        }

        public bool IsBefore(int pos)
        {
            return OneBasedEnd < pos;
        }

        public bool IsAfter(Exon segment)
        {
            return OneBasedStart > segment.OneBasedStart && OneBasedEnd > segment.OneBasedEnd && OneBasedStart > segment.OneBasedEnd;
        }

        public bool IsAfter(int pos)
        {
            return OneBasedStart > pos;
        }

        public bool Overlaps(Exon segment)
        {
            return !IsBefore(segment) && !IsAfter(segment);
        }

        public bool Overlaps(int pos)
        {
            return !IsBefore(pos) && !IsAfter(pos);
        }

        public bool Includes(Exon segment)
        {
            return OneBasedStart <= segment.OneBasedStart && OneBasedEnd >= segment.OneBasedEnd;
        }

        public bool Includes(long pos)
        {
            return OneBasedStart <= pos && OneBasedEnd >= pos;
        }

        #endregion Public Methods

    }
}
