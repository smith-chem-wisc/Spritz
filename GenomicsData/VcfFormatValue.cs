using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GenomicsData
{

    public class VcfFormatValue
    {
        public string format_type { get; set; }

        /// <summary>
        /// Number of comma-separated fields
        /// </summary>
        public int number { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public string type { get; set; }
    }
}
