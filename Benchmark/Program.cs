using Genomics;
using GenomicsData;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using UsefulProteomicsDatabases;

namespace Benchmark
{
    class Program
    {
        private static void Main(string[] args)
        {
            Dictionary<string, Chromosome> chroms = read_chroms();
            GeneModel gtf = GTF.ReadGenomeFeatures(@"D:\GRCh38.81\Homo_sapiens.GRCh38.81.gtf", chroms);
            //GeneModel gtf = GTF.ReadGenomeFeatures(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\sample.gtf", chroms);
            Console.WriteLine();
            Console.WriteLine(gtf.genes.Count + " genes");
            Console.WriteLine(gtf.genes.Sum(g => g.transcripts.Count) + " transcripts");
            Console.WriteLine(gtf.genes.Sum(g => g.transcripts.Sum(t => t.exons.Count)) + " all exons");
            Console.WriteLine(gtf.genes.Sum(g => g.exons.Count) + " unique exons");
            List<Sample> samples = VCF.ReadVCF(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\A549_sample.vcf", gtf.chromosomes);
            Console.WriteLine(samples[0].sequence_variants.Count + " variants");
            List<string> ensembl_seqs = ProteinDbLoader.LoadProteinFasta(@"D:\GRCh38.81\Homo_sapiens.GRCh38.pep.all.fa", false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression).Select(x => x.BaseSequence).ToList();
            Dictionary<string, int> formal_charges = Loaders.GetFormalChargesDictionary(Loaders.LoadPsiMod(Path.Combine(Environment.CurrentDirectory, "PSI-MOD.obo2.xml")));
            Loaders.LoadElements(Path.Combine(Environment.CurrentDirectory, "elements.dat"));
            Loaders.LoadUniprot(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges);
            List<string> trembl_seqs = ProteinDbLoader.LoadProteinXML(@"D:\human170717_uniprot-proteome%3AUP000005640.xml", false, PtmListLoader.ReadModsFromFile(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges), false, new List<string>(), out Dictionary<string, Modification> unknown).Select(x => x.BaseSequence).ToList();
            List<string> swiss_seqs = ProteinDbLoader.LoadProteinXML(@"D:\humanReviewed170717_uniprot-proteome%3AUP000005640+reviewed%3Ayes.xml", false, PtmListLoader.ReadModsFromFile(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges), false, new List<string>(), out Dictionary<string, Modification> unknown2).Select(x => x.BaseSequence).ToList();
            List<string> gtf_seqs = gtf.genes.SelectMany(g => g.transcripts).Where(x => x.start_codon_start >= 0 && x.stop_codon_start >= 0).Select(t => t.translate().BaseSequence).ToList();
            Console.WriteLine();
            Console.WriteLine("Comparison with Ensembl pep.all sequences");
            Console.WriteLine(ensembl_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            Console.WriteLine(ensembl_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in ensembl.pep.all");
            Console.WriteLine(gtf_seqs.Except(ensembl_seqs).Count().ToString() + " protein sequences only in gtf seqs");
            Console.WriteLine();
            Console.WriteLine("Comparison with Trembl sequences");
            Console.WriteLine(trembl_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            Console.WriteLine(trembl_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in Trembl");
            Console.WriteLine(gtf_seqs.Except(trembl_seqs).Count().ToString() + " protein sequences only in gtf seqs");
            Console.WriteLine();
            Console.WriteLine("Comparison with SwissProt sequences");
            Console.WriteLine(swiss_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            Console.WriteLine(swiss_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in Trembl");
            Console.WriteLine(gtf_seqs.Except(swiss_seqs).Count().ToString() + " protein sequences only in gtf seqs");
            Console.ReadKey();
        }

        private static Dictionary<string, Chromosome> read_chroms()
        {
            Dictionary<string, Chromosome> chroms = GenomeFasta.ReadGenomeFasta(@"D:\GRCh38.81\Homo_sapiens.GRCh38.dna.primary_assembly.fa");
            //Dictionary<string, Chromosome> chroms = GenomeFasta.ReadGenomeFasta(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\chr1_sample.fa");
            foreach (var x in chroms)
            {
                Console.WriteLine(x.Key + "\t" + x.Value.Length);
            }
            return chroms;
        }
    }
}
