using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Interval
        : IComparable<Interval>
    {
        #region Public Properties

        /// <summary>
        /// Chromosome name
        /// </summary>
        public string ChromosomeID { get; set; }

        /// <summary>
        /// Strand of the chromosome that this exon is on.
        /// </summary>
        public string Strand { get; set; }

        /// <summary>
        /// Start position.
        /// </summary>
        public long OneBasedStart { get; set; } = -1;

        /// <summary>
        /// End position.
        /// </summary>
        public long OneBasedEnd { get; set; } = -1;

        /// <summary>
        /// Parent of this interval
        /// </summary>
        public Interval Parent { get; set; }

        /// <summary>
        /// Type of interval
        /// </summary>
        public EffectType IntervalType { get; set; }

        /// <summary>
        /// Variants contained in this interval
        /// </summary>
        public List<Variant> Variants { get; set; } = new List<Variant>();

        #endregion Public Properties

        #region Constructor

        /// <summary>
        /// Constructs an interval from chromosome, strand, start, and end
        /// </summary>
        /// <param name="chromID"></param>
        /// <param name="strand"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        public Interval(Interval parent, string chromosomeID, string strand, long oneBasedStart, long oneBasedEnd)
        {
            this.ChromosomeID = chromosomeID;
            this.Strand = strand;
            this.OneBasedStart = oneBasedStart;
            this.OneBasedEnd = oneBasedEnd;
        }

        /// <summary>
        /// Constructs an interval
        /// </summary>
        /// <param name="interval"></param>
        public Interval(Interval interval) :
            this(interval.Parent, interval.ChromosomeID, interval.Strand, interval.OneBasedStart, interval.OneBasedEnd)
        {
        }

        public Interval()
        {
        }

        #endregion Constructor

        #region Variant Application Methods

        /// <summary>
        /// Do something with a variant within this interval
        /// </summary>
        /// <param name="variant"></param>
        public virtual Interval ApplyVariant(Variant variant)
        {
            Variants.Add(variant);

            Interval newInterval = null;
            switch (variant.VarType)
            {
                case Variant.VariantType.SNV:
                case Variant.VariantType.MNV:
                    // Variant does not change length. No effect when applying (all coordinates remain the same)
                    newInterval = this;
                    break;

                case Variant.VariantType.INS:
                    newInterval = ApplyIns(variant);
                    break;

                case Variant.VariantType.DEL:
                    newInterval = ApplyDel(variant);
                    break;

                case Variant.VariantType.DUP:
                    newInterval = ApplyDup(variant);
                    break;

                default:
                    // We are not ready for mixed changes
                    throw new ArgumentException("Variant type not supported: " + variant.VarType.ToString() + "\n\t" + variant);
            }

            // Always return a copy of the marker (if the variant is applied)
            if (newInterval == this)
            {
                return new Interval(this);
            }
            return newInterval;
        }

        /// <summary>
        /// Apply a Variant to a marker. Variant is a deletion
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        protected Interval ApplyDel(Variant variant)
        {
            Interval m = new Interval(this);

            if (variant.OneBasedEnd < m.OneBasedStart)
            {
                // Deletion before start: Adjust coordinates
                long lenChange = variant.LengthChange();
                m.OneBasedStart += lenChange;
                m.OneBasedEnd += lenChange;
            }
            else if (variant.Includes(m))
            {
                // Deletion completely includes this marker => The whole marker deleted
                return null;
            }
            else if (m.Includes(variant))
            {
                // This marker completely includes the deletion, but deletion does not include
                // marker. Marker is shortened (i.e. only 'end' coordinate needs to be updated)
                m.OneBasedEnd += variant.LengthChange();
            }
            else
            {
                // Variant is partially included in this marker.
                // This is treated as three different type of deletions:
                //		1- One after the marker
                //		2- One inside the marker
                //		3- One before the marker
                // Note that type 1 and 3 cannot exists at the same time, otherwise the
                // deletion would fully include the marker (previous case)

                // Part 1: Deletion after the marker
                if (m.OneBasedEnd < variant.OneBasedEnd)
                {
                    // Actually this does not affect the coordinates, so we don't care about this part
                }

                // Part 2: Deletion matching the marker (intersection)
                long istart = Math.Max(variant.OneBasedStart, m.OneBasedStart);
                long iend = Math.Min(variant.OneBasedEnd, m.OneBasedEnd);
                if (iend < istart) { throw new ArgumentOutOfRangeException("This should never happen!"); }// Sanity check
                m.OneBasedEnd -= (iend - istart + 1); // Update end coordinate

                // Part 3: Deletion before the marker
                if (variant.OneBasedStart < m.OneBasedEnd)
                {
                    // Update coordinates shifting the marker to the left
                    long delta = m.OneBasedStart - variant.OneBasedStart;
                    m.OneBasedStart -= delta;
                    m.OneBasedEnd -= delta;
                }
            }

            return m;
        }

        /// <summary>
        /// Apply a Variant to a marker. Variant is a duplication
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        protected Interval ApplyDup(Variant variant)
        {
            Interval m = new Interval(this);

            if (variant.OneBasedEnd < m.OneBasedStart)
            {
                // Duplication before marker start? => Adjust both coordinates
                long lenChange = variant.LengthChange();
                m.OneBasedStart += lenChange;
                m.OneBasedEnd += lenChange;
            }
            else if (variant.Includes(m))
            {
                // Duplication includes whole marker? => Adjust both coordinates
                long lenChange = variant.LengthChange();
                m.OneBasedStart += lenChange;
                m.OneBasedEnd += lenChange;
            }
            else if (m.Includes(variant))
            {
                // Duplication included in marker? => Adjust end coordinate
                m.OneBasedEnd += variant.LengthChange();
            }
            else if (variant.Intersects(m))
            {
                // Duplication includes part of marker? => Adjust end
                m.OneBasedEnd += variant.IntersectSize(m);
            }
            else
            {
                // Duplication after end, no effect on marker coordinates
            }

            return m;
        }

        /// <summary>
        /// Apply a Variant to a marker. Variant is an insertion
        /// </summary>
        /// <param name="variant"></param>
        /// <returns></returns>
        protected Interval ApplyIns(Variant variant)
        {
            Interval m = new Interval(this);

            if (variant.OneBasedStart < m.OneBasedStart)
            {
                // Insertion point before marker start? => Adjust both coordinates
                long lenChange = variant.LengthChange();
                m.OneBasedStart += lenChange;
                m.OneBasedEnd += lenChange;
            }
            else if (variant.OneBasedStart <= m.OneBasedEnd)
            {
                // Insertion point after start, but before end? => Adjust end coordinate
                m.OneBasedEnd += variant.LengthChange();
            }
            else
            {
                // Insertion point after end, no effect on marker coordinates
            }

            return m;
        }

        /// <summary>
        /// Calculate the effect of this variant
        /// </summary>
        /// <param name="variant"></param>
        /// <param name="variantEffects"></param>
        /// <returns></returns>
        public virtual bool VariantEffect(Variant variant, VariantEffects variantEffects)
        {
            if (!Intersects(variant)) { return false; }
            variantEffects.add(variant, this, IntervalType, "");
            return true;
        }

        /// <summary>
        /// Distance from the beginning/end of a list of intervals, until this SNP It count the number of bases in 'markers'
        /// </summary>
        /// <param name="markers"></param>
        /// <param name="fromEnd"></param>
        /// <returns></returns>
        public long DistanceBases(List<Interval> markers, bool fromEnd)
        {
            // Create a new list of sorted intervals
            List<Interval> markersSort = fromEnd ?
                new List<Interval>(markers).OrderBy(m => m.OneBasedEnd).ToList() :
                new List<Interval>(markers).OrderBy(m => m.OneBasedStart).ToList();

            // Calculate distance
            long len = 0, latest = -1;
            foreach (Interval m in markersSort)
            {
                // Initialize
                if (latest < 0)
                {
                    latest = fromEnd ? m.OneBasedEnd + 1 : m.OneBasedStart - 1;
                }

                if (fromEnd)
                {
                    if (Intersects(m)) { return len + (m.OneBasedEnd - OneBasedStart); }
                    else if (OneBasedStart > m.OneBasedEnd) { return len - 1 + (latest - OneBasedStart); }

                    latest = m.OneBasedStart;
                }
                else
                {
                    if (Intersects(m)) { return len + (OneBasedStart - m.OneBasedStart); }
                    else if (OneBasedStart < m.OneBasedStart) { return len - 1 + (OneBasedStart - latest); }

                    latest = m.OneBasedEnd;
                }

                len += m.Length();
            }

            return fromEnd ?
                len - 1 + (latest - OneBasedStart) :
                len - 1 + (OneBasedStart - latest);
        }

        #endregion Variant Application Methods

        #region Public Methods

        /// <summary>
        /// Compare by start and end
        /// </summary>
        /// <param name="i2"></param>
        /// <returns></returns>
        public virtual int CompareTo(Interval i2)
        {
            // Start
            if (OneBasedStart > i2.OneBasedStart) { return 1; }
            if (OneBasedStart < i2.OneBasedStart) { return -1; }

            // End
            if (OneBasedEnd > i2.OneBasedEnd) { return 1; }
            if (OneBasedEnd < i2.OneBasedEnd) { return -1; }

            return 0;
        }

        public Interval FindParent(Type type)
        {
            if (GetType().Equals(type))
            {
                return this;
            }
            if (Parent != null && Parent.GetType().Equals(type))
            {
                return Parent.FindParent(type);
            }
            return null;
        }

        /// <summary>
        /// Is this interval on the forward strand?
        /// </summary>
        /// <returns></returns>
        public bool IsStrandPlus()
        {
            return Strand == "+";
        }

        /// <summary>
        /// Is this interval on the reverse strand?
        /// </summary>
        /// <returns></returns>
        public bool IsStrandMinus()
        {
            return Strand != "+";
        }

        /// <summary>
        /// Determines whether this interval is before the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        public bool IsBefore(Interval segment)
        {
            return OneBasedStart < segment.OneBasedStart && OneBasedEnd < segment.OneBasedEnd && OneBasedEnd < segment.OneBasedStart;
        }

        /// <summary>
        /// Determines whether this interval is after the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        public bool IsAfter(Interval segment)
        {
            return OneBasedStart > segment.OneBasedStart && OneBasedEnd > segment.OneBasedEnd && OneBasedStart > segment.OneBasedEnd;
        }

        /// <summary>
        /// Determines whether this interval overlaps the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        public bool Intersects(Interval segment)
        {
            return !IsBefore(segment) && !IsAfter(segment);
        }

        /// <summary>
        /// Determines whether this interval includes the queried interval
        /// </summary>
        /// <param name="segment"></param>
        /// <returns></returns>
        public bool Includes(Interval segment)
        {
            return OneBasedStart <= segment.OneBasedStart && OneBasedEnd >= segment.OneBasedEnd;
        }

        /// <summary>
        /// Determines whether this interval is before the queried location
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public bool IsBefore(long pos)
        {
            return OneBasedEnd < pos;
        }

        /// <summary>
        /// Determines whether this interval is after the queried location
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public bool IsAfter(long pos)
        {
            return OneBasedStart > pos;
        }

        /// <summary>
        /// Determines whether this interval contains the queried location
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public bool Intersects(long pos)
        {
            return !IsBefore(pos) && !IsAfter(pos);
        }

        /// <summary>
        /// Determines whether this interval contains the queried location
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public bool Includes(long pos)
        {
            return OneBasedStart <= pos && OneBasedEnd >= pos;
        }

        public long Length()
        {
            return OneBasedEnd - OneBasedStart + 1;
        }

        public Interval Intersect(Interval interval)
        {
            if (interval.ChromosomeID != ChromosomeID) { return null; }
            long istart = Math.Max(OneBasedStart, interval.OneBasedStart);
            long iend = Math.Min(OneBasedEnd, interval.OneBasedEnd);
            if (iend < istart) { return null; }
            return new Interval(this);
        }

        public long IntersectSize(Interval other)
        {
            long start = Math.Max(this.OneBasedStart, other.OneBasedStart);
            long end = Math.Min(this.OneBasedEnd, other.OneBasedEnd);
            if (end < start) { return 0; }
            return end - start + 1;
        }

        #endregion Public Methods

        #region Public Static Method

        public static long GetMedian(IEnumerable<Interval> intervals)
        {
            // Add all start and end coordinates
            long i = 0;
            long[] points = new long[2 * intervals.Count()];
            foreach (Interval interval in intervals)
            {
                points[i++] = interval.OneBasedStart;
                points[i++] = interval.OneBasedEnd;
            }

            Array.Sort(points);
            int middle = points.Length / 2;
            return points[middle];
        }

        #endregion Public Static Method
    }
}