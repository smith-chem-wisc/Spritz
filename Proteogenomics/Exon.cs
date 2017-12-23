using Bio;
using Bio.VCF;
using Bio.Extensions;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System;

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

        public List<Exon> GetExonSequences(int maxCombos, bool getVariantSequences = false, double homozygousThreshold = 1, bool phased = false, int fakePhasingRange = 101)
        {
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

            // This could be VERY inefficient, e.g. even with just 24 adjacent variants
            //int maxCount = (int)Math.Log(Int32.MaxValue, 2);
            List<List<Variant>> possibleHaplotypes = new List<List<Variant>> { new List<Variant>(homozygous) };
            possibleHaplotypes.AddRange(
                Enumerable.Range(1, heterozygous.Count).SelectMany(k =>
                    ExtensionMethods.Combinations(heterozygous, k)
                        .Select(heteroAlleles => new List<Variant>(homozygous).Concat(heteroAlleles).ToList()))
                    .ToList());

            // fake phasing (of variants within 100 bp of each other) to eliminate combinitorial explosion
            // this section eliminates instances where variants that are very close to each other are emitted as separate haplotypes
            if (!phased)
            {
                int range = fakePhasingRange;
                List<List<Variant>> fakePhased = null;
                while (fakePhased == null || fakePhased.Count > Math.Max(2, maxCombos))
                {
                    fakePhased = new List<List<Variant>>();
                    List<List<Variant>> phasedLast = possibleHaplotypes.OrderBy(x => x.Count).ToList();
                    List<int> variantSites = phasedLast.Last().Select(v => v.OneBasedStart).ToList();
                    foreach (List<Variant> hap in phasedLast)
                    {
                        List<int> hapSites = hap.Select(v => v.OneBasedStart).ToList();
                        bool containsVariantsToFakePhase = variantSites.Any(vs => !hapSites.Contains(vs) && hap.Any(v => Math.Abs(v.OneBasedStart - vs) < range));
                        if (!containsVariantsToFakePhase)
                            fakePhased.Add(hap);
                    }
                    possibleHaplotypes = fakePhased;
                    range *= 2;
                }
            }

            // get internal byte array for sequence copies
            // using byte arrays for intermediate data structures greatly improves performance because of copying and storage efficiency over strings
            byte[] sequenceBytes = new byte[this.Sequence.Count];
            (this.Sequence as Sequence).CopyTo(sequenceBytes, 0, this.Sequence.Count);

            List<Exon> sequences = new List<Exon>();
            foreach (List<Variant> haplotype in possibleHaplotypes)
            {
                long tmpOneBasedStart = OneBasedStart;
                byte[] newSequence = new byte[this.Sequence.Count + haplotype.Sum(v => v.AlternateAllele.Length - 1)];
                haplotype.OrderBy(v => v.OneBasedStart);

                int subseqStart;
                int subseqCount;
                foreach (Variant var in haplotype)
                {
                    if (var.OneBasedStart < tmpOneBasedStart)
                        continue; // variant collides with previous variant

                    subseqStart = (int)(tmpOneBasedStart - this.OneBasedStart);
                    subseqCount = (int)(var.OneBasedStart - tmpOneBasedStart);
                    if (subseqStart + subseqCount - 1 > this.Sequence.Count) // allele extends beyond exon, so mind the exon boundary
                    {
                        subseqCount = (int)(this.Sequence.Count - subseqStart);
                        long subAlleleCount = this.OneBasedEnd - var.OneBasedStart;
                        Array.Copy(sequenceBytes, subseqStart, newSequence, subseqStart, subseqCount);
                        Array.Copy(Encoding.UTF8.GetBytes(var.AlternateAllele), 0, newSequence, subseqStart + subseqCount, subAlleleCount);
                        long newLength = newSequence.Length + this.OneBasedEnd - var.AlternateAllele.Length;
                        byte[] tmpNewSequence = new byte[newLength];
                        Array.Copy(newSequence, tmpNewSequence, newLength);
                        newSequence = tmpNewSequence;
                    }
                    else
                    {
                        Array.Copy(sequenceBytes, subseqStart, newSequence, subseqStart, subseqCount);
                        Array.Copy(Encoding.UTF8.GetBytes(var.AlternateAllele), 0, newSequence, subseqStart + subseqCount, var.AlternateAllele.Length);
                    }
                    tmpOneBasedStart = var.OneBasedStart + var.ReferenceAllele.Length;
                }

                // Append the remaining sequence
                subseqStart = (int)(tmpOneBasedStart - this.OneBasedStart);
                subseqCount = (int)(OneBasedEnd - tmpOneBasedStart + 1);
                Array.Copy(sequenceBytes, subseqStart, newSequence, subseqStart, subseqCount);
                Sequence seq = new Sequence(this.Sequence.Alphabet, newSequence);
                Exon x = new Exon(seq, this.OneBasedStart, this.OneBasedEnd, this.ChromID, this.Metadata);
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
