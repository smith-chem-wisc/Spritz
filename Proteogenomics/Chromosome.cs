using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Bio;

namespace Proteogenomics
{
    public class Chromosome
        : Interval
    {

        public ISequence Sequence { get; set; }

        public Chromosome(string chromID, long oneBasedStart, long oneBasedEnd)
            : base(chromID, "+", oneBasedStart, oneBasedEnd)
        {

        }
    }
}
