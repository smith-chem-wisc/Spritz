using Bio;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Interval
        : IComparable<Interval>
    {
        /// <summary>
        /// Constructs an interval from chromosome, strand, start, and end
        /// </summary>
        /// <param name="chromID"></param>
        /// <param name="strand"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        public Interval(Interval parent, string chromosomeID, string source, string strand, long oneBasedStart, long oneBasedEnd)
        {
            Parent = parent;
            ChromosomeID = chromosomeID;
            Source = source;
            Strand = strand;
            OneBasedStart = oneBasedStart;
            OneBasedEnd = oneBasedEnd;
        }

        /// <summary>
        /// Constructs an interval
        /// </summary>
        /// <param name="interval"></param>
        public Interval(Interval interval) :
            this(interval.Parent, interval.ChromosomeID, interval.Source, interval.Strand, interval.OneBasedStart, interval.OneBasedEnd)
        {
        }

        public Interval()
        {
        }

        /// <summary>
        /// Chromosome name
        /// </summary>
        public string ChromosomeID { get; set; }

        /// <summary>
        /// Source name
        /// </summary>
        public string Source { get; set; }

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
        /// Feature name used for writing GTF files
        /// </summary>
        public virtual string FeatureType { get; } = "interval";

        /// <summary>
        /// Get the median of a set of intervals
        /// </summary>
        /// <param name="intervals"></param>
        /// <returns></returns>
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
                new List<Interval>(markers).OrderByDescending(m => m.OneBasedEnd).ToList() :
                new List<Interval>(markers).OrderBy(m => m.OneBasedStart).ToList();

            // Calculate distance
            long len = 0;
            long latest = -1;
            foreach (Interval m in markersSort)
            {
                // Initialize
                if (latest < 0)
                {
                    latest = fromEnd ?
                        m.OneBasedEnd + 1 :
                        m.OneBasedStart - 1;
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

        public List<Interval> Minus(Interval interval)
        {
            List<Interval> intervals = new List<Interval>();
            if (Intersects(interval))
            {
                if (interval.Includes(this))
                {
                    // 'this' is included in 'interval' => Nothing left
                }
                else if (interval.OneBasedStart <= OneBasedStart && interval.OneBasedEnd < OneBasedEnd)
                {
                    // 'interval' overlaps left part of 'this' => Include right part of 'this'
                    intervals.Add(new Interval(Parent, ChromosomeID, Source, Strand, interval.OneBasedEnd + 1, OneBasedEnd));
                }
                else if (OneBasedStart < interval.OneBasedStart && OneBasedEnd <= interval.OneBasedEnd)
                {
                    // 'interval' overlaps right part of 'this' => Include left part of 'this'
                    intervals.Add(new Interval(Parent, ChromosomeID, Source, Strand, OneBasedStart, interval.OneBasedStart - 1));
                }
                else if (OneBasedStart < interval.OneBasedStart && interval.OneBasedEnd < OneBasedEnd)
                {
                    // 'interval' overlaps middle of 'this' => Include left and right part of 'this'
                    intervals.Add(new Interval(Parent, ChromosomeID, Source, Strand, OneBasedStart, interval.OneBasedStart - 1));
                    intervals.Add(new Interval(Parent, ChromosomeID, Source, Strand, interval.OneBasedEnd + 1, OneBasedEnd));
                }
                else
                {
                    throw new ArgumentException("Interval intersection not analysed. This should nbever happen!");
                }
            }
            else
            {
                intervals.Add(this); // No intersection => Just add 'this' interval
            }

            return intervals;
        }

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
            if (Parent != null && Parent as Interval != null)
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
            return segment.OneBasedEnd >= OneBasedStart && segment.OneBasedStart <= OneBasedEnd;
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
            return OneBasedStart <= pos && pos <= OneBasedEnd;
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

        public MetadataListItem<List<string>> GetGtfFeatureMetadata()
        {
            var feature = new MetadataListItem<List<string>>(FeatureType, GetGtfAttributes());
            feature.SubItems["source"] = new List<string> { Source.ToString() };
            feature.SubItems["start"] = new List<string> { OneBasedStart.ToString() };
            feature.SubItems["end"] = new List<string> { OneBasedEnd.ToString() };
            if (Strand != ".") { feature.SubItems["strand"] = new List<string> { Strand.ToString() }; } // might take in features without strand later on
            return feature;
        }

        /// <summary>
        /// Gets the attributes string as free text (default is empty string)
        /// </summary>
        /// <returns></returns>
        public virtual string GetGtfAttributes()
        {
            return "";
        }
    }
}