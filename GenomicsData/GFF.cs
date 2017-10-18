using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using Bio.IO.Gff;

namespace GenomicsData
{
    public class GFF
    {
        public GFF(string filename)
        {
            GffParser gff = new GffParser();
            var x = gff.Parse(new FileStream(filename, FileMode.Open));
            (x.FirstOrDefault() as Bio.Sequence).
        }
    }
}
