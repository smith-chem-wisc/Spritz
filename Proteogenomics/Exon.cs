using Bio;
using Bio.VCF;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class Exon :
        Interval
    {

        #region Public Properties

        /// <summary>
        /// Sequence using BioDotNet interface
        /// </summary>
        public ISequence Sequence { get; set; }

        /// <summary>
        /// Used for both objects straight from VCF and for specific Variants, which extend VariantContext
        /// </summary>
        public List<VariantContext> Variants { get; set; } = new List<VariantContext>();

        #endregion Public Properties

        #region Public Constructor

        public Exon(ISequence Sequence, long oneBasedStart, long oneBasedEnd, string chromID, string strand)
            : base(chromID, strand, oneBasedStart, oneBasedEnd)
        {
            this.Sequence = Sequence;
        }

        #endregion Public Constructor

        #region Private Constructor

        private Exon(Exon exon)
        {
            this.Sequence = exon.Sequence;
            this.OneBasedStart = exon.OneBasedStart;
            this.OneBasedEnd = exon.OneBasedEnd;
            this.ChromosomeID = exon.ChromosomeID;
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
                    ProteogenomicsUtility.Combinations(heterozygous, k)
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
                Exon x = new Exon(seq, this.OneBasedStart, this.OneBasedEnd, this.ChromosomeID, this.Metadata);
                x.Variants = haplotype.OfType<VariantContext>().ToList();
                sequences.Add(x);
            }
            return sequences;
        }

        #endregion Public Methods

    }
}
