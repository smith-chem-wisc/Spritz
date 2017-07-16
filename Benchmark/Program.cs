using System.Collections.Generic;
using System.Linq;
using Genomics;
using GenomicsData;
using System;

namespace Benchmark
{
    class Program
    {
        private static void Main(string[] args)
        {
            Dictionary<string, Chromosome> chroms = read_chroms();
            //GeneModel gtf = GTF.ReadGenomeFeatures(@"D:\GRCh38.81\Homo_sapiens.GRCh38.81.gtf", chroms);
            GeneModel gtf = GTF.ReadGenomeFeatures(@"D:\GRCh38.81\sample.gtf", chroms);
            Console.WriteLine();
            Console.WriteLine(gtf.genes.Count + " genes");
            Console.WriteLine(gtf.genes.Sum(g => g.transcripts.Count) + " transcripts");
            Console.WriteLine(gtf.genes.Sum(g => g.transcripts.Sum(t => t.exons.Count)) + " all exons");
            Console.WriteLine(gtf.genes.Sum(g => g.exons.Count) + " unique exons");
            List<Sample> samples = VCF.ReadVCF(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\A549_sample.vcf", gtf.chromosomes);
            Console.WriteLine(samples[0].sequence_variants.Count + " variants");
            foreach (Transcript t in gtf.genes.SelectMany(g => g.transcripts).Where(x => x.start_codon_start >= 0 && x.stop_codon_start >=0))
            {
                Console.WriteLine(t.translate().BaseSequence);
            }
            Console.ReadKey();
        }

        private static Dictionary<string, Chromosome> read_chroms()
        {
            //Dictionary<string, Chromosome> chroms = GenomeFasta.ReadGenomeFasta(@"D:\GRCh38.81\Homo_sapiens_GRCh38.dna.primary.assembly.canonPlus.fa");
            Dictionary<string, Chromosome> chroms = GenomeFasta.ReadGenomeFasta(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\chr1_sample.fa");
            foreach (var x in chroms)
            {
                Console.WriteLine(x.Key + "\t" + x.Value.Length);
            }
            return chroms;
        }
    }
}
