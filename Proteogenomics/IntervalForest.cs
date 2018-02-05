using System.Collections.Generic;

namespace Proteogenomics
{
    /// <summary>
    /// A set of interval trees, e.g. one per chromosome, one per transcript, etc.
    /// </summary>
    public class IntervalForest 
    {

        #region Public Properties

        public string Name { get; set; }

        public Dictionary<string, IntervalTree> Forest { get; set; } = new Dictionary<string, IntervalTree>();

        #endregion Public Properties

        #region Constructors

        public IntervalForest()
        {
        }

        public IntervalForest(IEnumerable<Interval> intervals)
        {
            Add(intervals);
        }

        #endregion Constructors

        #region Public Methods

        public void Add(IEnumerable<Interval> intervals)
        {
            foreach (Interval i in intervals)
            {
                Add(i);
            }
        }

        public void Add(Interval interval)
        {
            if (interval == null) return;
            if (Forest.TryGetValue(interval.ChromID, out IntervalTree tree)) tree.Add(interval);
            else Forest.Add(interval.ChromID, new IntervalTree(new List<Interval> { interval }));
        }

        public void Build()
        {
            foreach (IntervalTree it in Forest.Values)
            {
                it.Build();
            }
        }

        #endregion Public Methods
    }
}
