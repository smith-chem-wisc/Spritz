using System.Collections.Generic;
using System.Linq;

namespace Genomics
{
    public class ChromosomeSegment
    {
        #region Public Properties

        public string ID { get; set; }
        public Chromosome Chrom { get; set; }
        public string Strand { get; set; }
        public int OneBasedStart { get; set; }
        public int OneBasedEnd { get; set; }
        public int Length { get { return OneBasedEnd - OneBasedStart; } }

        #endregion Public Properties

        #region Public Constructor

        public ChromosomeSegment(string id, Chromosome chrom, string strand, int ZeroBasedStart, int ZeroBasedEnd)
        {
            this.ID = id;
            this.Chrom = chrom;
            this.Strand = strand;
            this.OneBasedStart = ZeroBasedStart;
            this.OneBasedEnd = ZeroBasedEnd;
        }

        #endregion Public Constructor

        #region Public Comparison Methods

        public bool is_before(ChromosomeSegment segment)
        {
            return OneBasedStart < segment.OneBasedStart && OneBasedEnd < segment.OneBasedEnd && OneBasedEnd < segment.OneBasedStart;
        }

        public bool is_before(int pos)
        {
            return OneBasedEnd < pos;
        }

        public bool is_after(ChromosomeSegment segment)
        {
            return OneBasedStart > segment.OneBasedStart && OneBasedEnd > segment.OneBasedEnd && OneBasedStart > segment.OneBasedEnd;
        }

        public bool is_after(int pos)
        {
            return OneBasedStart > pos;
        }

        public bool overlaps(ChromosomeSegment segment)
        {
            return !is_before(segment) && !is_after(segment);
        }

        public bool overlaps(int pos)
        {
            return !is_before(pos) && !is_after(pos);
        }

        public bool includes(ChromosomeSegment segment)
        {
            return OneBasedStart <= segment.OneBasedStart && OneBasedEnd >= segment.OneBasedEnd;
        }

        public bool includes(int pos)
        {
            return Enumerable.Range(OneBasedStart, OneBasedEnd - OneBasedStart + 1).Contains(pos);
        }

        public bool equals(ChromosomeSegment segment)
        {
            return Chrom == segment.Chrom && OneBasedStart == segment.OneBasedStart && OneBasedEnd == segment.OneBasedEnd;
        }

        public bool equals(int pos)
        {
            return OneBasedStart == OneBasedEnd && OneBasedStart == pos;
        }

        #endregion Public Comparison Methods

        #region Public Methods

        //public ChromosomeSegment combine_overlapping_sequences(IEnumerable<ChromosomeSegment> other_seqs, out IEnumerable<ChromosomeSegment> combined_with_this)
        //{

        //}

        public override bool Equals(object obj)
        {
            ChromosomeSegment seg = obj as ChromosomeSegment;
            return
                seg != null &&
                seg.Chrom.Name == Chrom.Name &&
                seg.Strand == Strand &&
                seg.OneBasedStart == OneBasedStart &&
                seg.OneBasedEnd == OneBasedEnd;
        }

        public override int GetHashCode()
        {
            return Chrom.GetHashCode() ^ Strand.GetHashCode() ^ OneBasedStart ^ OneBasedEnd;
        }

        #endregion Public Methods
    }
}
