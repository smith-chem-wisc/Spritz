using Bio;
using Bio.VCF;
using Bio.Extensions;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Genomics
{
    public class Exon
    {

        #region Private Properties

        private MetadataListItem<List<string>> metadata { get; set; }

        #endregion Private Properties

        #region Public Properties

        public ISequence Sequence { get; set; }
        public string ChromID { get; set; }
        public string Strand { get; set; }
        public long OneBasedStart { get; set; } = -1;
        public long OneBasedEnd { get; set; } = -1;
        public List<VariantContext> variants { get; set; } = new List<VariantContext>();

        #endregion Public Properties

        #region Public Constructor

        public Exon(ISequence Sequence, long start, long stop, string chrom_id, MetadataListItem<List<string>> metadata)
        {
            this.Sequence = Sequence;
            this.OneBasedStart = start;
            this.OneBasedEnd = stop;
            this.ChromID = chrom_id;
            this.Strand = metadata.SubItems["strand"][0];
            this.metadata = metadata;
        }

        #endregion Public Constructor

        #region Public Methods

        public IEnumerable<string> get_sequences(double homozygous_threshold, bool get_variant_sequences)
        {
            string seq = SequenceExtensions.ConvertToString(this.Sequence);
            if (variants.Count == 0 || !get_variant_sequences)
                return new List<string> { seq };

            List<VariantLite> homozygous = new List<VariantLite>();
            List<VariantLite> heterozygous = new List<VariantLite>();
            foreach (VariantLite var in variants.SelectMany(v => VariantLite.ParseVariantContext(v)))
            {
                // todo: allow depth filter here using "DP" column, especially for indels
                if (var.AlleleFrequency >= homozygous_threshold) homozygous.Add(var); // homozygous if the allele frequency is about 1
                else heterozygous.Add(var);
            }

            List<List<VariantLite>> possible_haplotypes = Enumerable.Range(1, heterozygous.Count).SelectMany(k =>
                ExtensionMethods.Combinations(heterozygous, k)
                    .Select(hetero_alleles => new List<VariantLite>(homozygous).Concat(hetero_alleles).ToList()))
                .ToList();

            List<string> sequences = new List<string>();
            foreach (List<VariantLite> haplotype in possible_haplotypes)
            {
                long tmpOneBasedStart = OneBasedStart;
                StringBuilder builder = new StringBuilder();
                haplotype.OrderBy(v => v.OneBasedStart);
                foreach (VariantLite var in haplotype)
                {
                    builder.Append(seq.Substring((int)(tmpOneBasedStart - OneBasedStart), (int)(var.OneBasedStart - this.OneBasedStart)));
                    builder.Append(var.AlternateAllele);
                    tmpOneBasedStart = var.OneBasedStart + var.ReferenceAllele.Length;
                }
                sequences.Add(builder.ToString());
            }
            return sequences;
        }

        //public override int GetHashCode()
        //{
        //    return metadata.Key ^ metadata.SubItems["features"];
        //}

        public bool is_before(Exon segment)
        {
            return OneBasedStart < segment.OneBasedStart && OneBasedEnd < segment.OneBasedEnd && OneBasedEnd < segment.OneBasedStart;
        }

        public bool is_before(int pos)
        {
            return OneBasedEnd < pos;
        }

        public bool is_after(Exon segment)
        {
            return OneBasedStart > segment.OneBasedStart && OneBasedEnd > segment.OneBasedEnd && OneBasedStart > segment.OneBasedEnd;
        }

        public bool is_after(int pos)
        {
            return OneBasedStart > pos;
        }

        public bool overlaps(Exon segment)
        {
            return !is_before(segment) && !is_after(segment);
        }

        public bool overlaps(int pos)
        {
            return !is_before(pos) && !is_after(pos);
        }

        public bool includes(Exon segment)
        {
            return OneBasedStart <= segment.OneBasedStart && OneBasedEnd >= segment.OneBasedEnd;
        }

        public bool includes(long pos)
        {
            return OneBasedStart <= pos && OneBasedEnd >= pos;
        }

        #endregion Public Methods

    }
}
