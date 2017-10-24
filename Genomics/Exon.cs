using Bio;
using System;
using System.Collections.Generic;
using System.Linq;
using Bio.VCF;

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
        List<VariantContext> variants { get; set; } = new List<VariantContext>();

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
