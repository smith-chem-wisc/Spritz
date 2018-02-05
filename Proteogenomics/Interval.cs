using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class Interval
    {

        #region Public Properties

        /// <summary>
        /// Chromosome name
        /// </summary>
        public string ChromID { get; set; }

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

        #endregion Public Properties

        #region Constructor

        /// <summary>
        /// Constructs an interval from chromosome, strand, start, and end
        /// </summary>
        /// <param name="chromID"></param>
        /// <param name="strand"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        public Interval(string chromID, string strand, long oneBasedStart, long oneBasedEnd)
        {
            this.ChromID = chromID;
            this.Strand = strand;
            this.OneBasedStart = oneBasedStart;
            this.OneBasedEnd = oneBasedEnd;
        }

        /// <summary>
        /// Constructs an interval
        /// </summary>
        /// <param name="interval"></param>
        public Interval(Interval interval) :
            this(interval.ChromID, interval.Strand, interval.OneBasedStart, interval.OneBasedEnd)
        {

        }

        #endregion Constructor

        #region Public Methods

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
        public bool Overlaps(Interval segment)
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
        public bool Overlaps(long pos)
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

        #endregion Public Methods

        #region Public Static Method

        public static long GetMedian(List<Interval> intervals)
        {
            // Add all start and end coordinates
            long i = 0;
            long[] points = new long[2 * intervals.Count];
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
