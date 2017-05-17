using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Genomics
{
    public class ChromosomeSegment
    {
        public string name { get; set; }
        public string id { get; set; }
        public string biotype { get; set; }
        public Chromosome chrom { get; set; }
        public string strand { get; set; }
        public int start { get; set; }
        public int end { get; set; }
        public int length { get { return this.end - this.start; } }

        public ChromosomeSegment(string id, Chromosome chrom, string strand, int start, int end, string name, string biotype)
        {
            this.name = name;
            this.id = id;
            this.biotype = biotype;
            this.chrom = chrom;
            this.strand = strand;
            this.start = start;
            this.end = end;
        }

        public bool is_before(ChromosomeSegment segment)
        {
            return this.start < segment.start && this.end < segment.end && this.end < segment.start;
        }
        public bool is_before(int pos)
        {
            return this.end < pos;
        }

        public bool is_after(ChromosomeSegment segment)
        {
            return this.start > segment.start && this.end > segment.end && this.start > segment.end;
        }
        public bool is_after(int pos)
        {
            return this.start > pos;
        }

        public bool overlaps(ChromosomeSegment segment)
        {
            return !this.is_before(segment) && !this.is_after(segment);
        }
        public bool overlaps(int pos)
        {
            return !this.is_before(pos) && !this.is_after(pos);
        }

        public bool includes(ChromosomeSegment segment)
        {
            return this.start <= segment.start && this.end >= segment.end;
        }
        public bool includes(int pos)
        {
            return Enumerable.Range(this.start, this.end - this.start + 1).Contains(pos);
        }

        public bool equals(ChromosomeSegment segment)
        {
            return this.chrom == segment.chrom && this.start == segment.start && this.end == segment.end;
        }
        public bool equals(int pos)
        {
            return this.start == this.end && this.start == pos;
        }
    }
}
