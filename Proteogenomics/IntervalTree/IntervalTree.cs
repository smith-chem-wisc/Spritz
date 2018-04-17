using System.Collections.Generic;

namespace Proteogenomics
{
    public class IntervalTree
    {
        public IntervalNode Head { get; set; } = new IntervalNode();

        public List<Interval> Intervals { get; set; } = new List<Interval>();

        public bool Synced { get; set; }

        /// <summary>
        /// Instantiate a new interval tree with no intervals
        /// </summary>
        public IntervalTree()
        {
        }

        /// <summary>
        /// Instantiate a new interval tree with a list of intervals
        /// </summary>
        /// <param name="intervals"></param>
        public IntervalTree(IEnumerable<Interval> intervals)
        {
            Head = new IntervalNode(intervals);
            Intervals = new List<Interval>(intervals);
            Synced = true;
        }

        public void Add(Interval interval)
        {
            Intervals.Add(interval);
            Synced = false;
        }

        public void Add(IEnumerable<Interval> intervals)
        {
            Intervals.AddRange(intervals);
            Synced = false;
        }

        public void Build()
        {
            if (!Synced)
            {
                lock (this)
                {
                    Head = new IntervalNode(Intervals);
                    Synced = true;
                }
            }
        }

        public void Build(IEnumerable<Interval> intervals)
        {
            Add(intervals);
            Build();
        }

        public List<Interval> Query(Interval interval)
        {
            if (!Synced) { Build(); }
            return Head.Query(interval);
        }

        public List<Interval> Stab(long point)
        {
            if (!Synced) { Build(); }
            return Head.Stab(point);
        }
    }
}