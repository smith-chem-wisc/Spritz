using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Genomics;
using GenomicsData;

namespace Benchmark
{
    class Program
    {
        private static void Main(string[] args)
        {
            VCF vcf = new VCF(@"E:\lelantos\linux_docs\TenCellLines\VCFs\A549.vcf");
            Console.WriteLine(vcf.samples[0].sequence_variants.Count);
            Console.ReadKey();
        }
    }
}
