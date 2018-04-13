using Bio;
using Bio.IO.FastA;
using Bio.IO.Gff;
using Bio.VCF;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using Proteogenomics;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Benchmark
{
    class Program
    {
        private static void Main(string[] args)
        {
            //Genome genome = new Genome(@"D:\GRCh38.81\Homo_sapiens.GRCh38.dna.primary_assembly.fa");
            //Console.WriteLine(genome.chroms.Count + " chroms");

            //GeneModel geneModel = new GeneModel(genome, @"D:\GRCh38.81\Homo_sapiens.GRCh38.81.gtf");
            //Console.WriteLine(geneModel.genes.Count + " genes");
            //Console.WriteLine(geneModel.genes.Sum(g => g.transcripts.Count) + " transcripts");
            //Console.WriteLine(geneModel.genes.Sum(g => g.exons.Count) + " exons");

            //List<Protein> proteins = geneModel.genes.SelectMany(g => g.translate()).ToList();
            //Console.WriteLine(proteins.Count + " proteins");

            VCFParser vcf = new VCFParser(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\A549_sample.vcf");
            List<VariantContext> variants = vcf.Select(x => x).ToList();
            Console.WriteLine(variants.Count + " alleles");

            Genome genome = new Genome(@"D:\GRCh37_canon\GRCh37_canon.fa");
            Console.WriteLine(genome.Chromosomes.Count + " chroms");

            List<Protein> ensembl_seqs = ProteinDbLoader.LoadProteinFasta(@"D:\GRCh37.87\Homo_sapiens.GRCh37.pep.all.fa", true, DecoyType.None, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression, ProteinDbLoader.uniprot_organism_expression, out List<string> errors).ToList();
            HashSet<string> incompletes = new HashSet<string>(ensembl_seqs.Where(p => 
                p.BaseSequence.StartsWith("X") || p.BaseSequence.EndsWith("X") || p.BaseSequence.Contains("*")).Select(p => p.Accession.Split(' ')[0].Split('.')[0]));

            GeneModel ensemblModel = new GeneModel(genome, @"D:\GRCh37.87\Homo_sapiens.GRCh37.87.gff3");
            Console.WriteLine(ensemblModel.Genes.Count + " genes");
            Console.WriteLine(ensemblModel.Genes.Sum(g => g.Transcripts.Count) + " transcripts");
            Console.WriteLine(ensemblModel.Genes.Sum(g => g.Transcripts.Sum(t => t.Exons.Count)) + " exons");
            List<Protein> ensembl_proteins = ensemblModel.Genes.SelectMany(g => g.Translate(true)).Where(p => !incompletes.Contains(p.Accession)).ToList();
            Console.WriteLine(ensembl_proteins.Count + " proteins");

            GeneModel pacBioModel = new GeneModel(genome, @"D:\GRCh37.87\IsoSeq_MCF7_2015edition_polished.unimapped.gff");
            Console.WriteLine(pacBioModel.Genes.Count + " genes2");
            Console.WriteLine(pacBioModel.Genes.Sum(g => g.Transcripts.Count) + " transcripts2");
            Console.WriteLine(pacBioModel.Genes.Sum(g => g.Transcripts.Sum(t => t.Exons.Count)) + " exons2");
            List<Protein> pacbio_proteins = pacBioModel.Genes.SelectMany(g => g.TranslateUsingAnnotatedStartCodons(ensemblModel, null, 7)).ToList();
            Dictionary<string, List<Protein>> seq_unique = new Dictionary<string, List<Protein>>();
            foreach (Protein p in pacbio_proteins)
            {
                seq_unique.TryGetValue(p.BaseSequence, out List<Protein> ps);
                if (ps == null) seq_unique.Add(p.BaseSequence, new List<Protein> { p });
                else ps.Add(p);
            }
            List<Protein> uniq_pacbio_proteins = seq_unique.Select(kv => kv.Value.First()).OrderBy(p => p.Accession).ToList();
            Console.WriteLine(uniq_pacbio_proteins.Count + " proteins");

            Console.WriteLine(ensembl_proteins.Select(p => p.BaseSequence).Intersect(uniq_pacbio_proteins.Select(p => p.BaseSequence)).Count().ToString() + " protein sequences in common");
            Console.WriteLine(ensembl_proteins.Select(p => p.BaseSequence).Except(uniq_pacbio_proteins.Select(p => p.BaseSequence)).Count().ToString() + " protein sequences only in ensembl");
            Console.WriteLine(uniq_pacbio_proteins.Select(p => p.BaseSequence).Except(ensembl_proteins.Select(p => p.BaseSequence)).Count().ToString() + " protein sequences only in pacbio");

            //GeneModel gtf = GTF.ReadGenomeFeatures(@"D:\GRCh38.81\Homo_sapiens.GRCh38.81.gtf", chroms);
            ////GeneModel gtf = GTF.ReadGenomeFeatures(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\sample.gtf", chroms);
            //Console.WriteLine();
            //Console.WriteLine(genes.Sum(g => g.transcripts.Count) + " transcripts");
            //Console.WriteLine(enes.Sum(g => g.transcripts.Sum(t => t.exons.Count)) + " all exons");
            //Console.WriteLine(gtf.genes.Sum(g => g.exons.Count) + " unique exons");
            //List<Sample> samples = VCF.ReadVCF(@"C:\Users\antho\Documents\GitHub\ProteoformDatabaseEngine\Test\A549_sample.vcf", gtf.chromosomes);
            //Console.WriteLine(samples[0].sequence_variants.Count + " variants");
            //List<string> ensembl_seqs = ProteinDbLoader.LoadProteinFasta(@"D:\GRCh38.81\Homo_sapiens.GRCh38.pep.all.fa", false, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression).Select(x => x.BaseSequence).ToList();
            //Dictionary<string, int> formal_charges = Loaders.GetFormalChargesDictionary(Loaders.LoadPsiMod(Path.Combine(Environment.CurrentDirectory, "PSI-MOD.obo2.xml")));
            //Loaders.LoadElements(Path.Combine(Environment.CurrentDirectory, "elements.dat"));
            //Loaders.LoadUniprot(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges);
            //List<string> trembl_seqs = ProteinDbLoader.LoadProteinXML(@"D:\human170717_uniprot-proteome%3AUP000005640.xml", true, DecoyType.None, PtmListLoader.ReadModsFromFile(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges), false, new List<string>(), out Dictionary<string, Modification> unknown).Select(x => x.BaseSequence).ToList();
            //List<string> swiss_seqs = ProteinDbLoader.LoadProteinXML(@"D:\humanReviewed170717_uniprot-proteome%3AUP000005640+reviewed%3Ayes.xml", true, DecoyType.None, PtmListLoader.ReadModsFromFile(Path.Combine(Environment.CurrentDirectory, "ptmlist.txt"), formal_charges), false, new List<string>(), out Dictionary<string, Modification> unknown2).Select(x => x.BaseSequence).ToList();
            //Console.WriteLine();
            //Console.WriteLine("Comparison with Ensembl pep.all sequences");
            //Console.WriteLine(ensembl_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            //Console.WriteLine(ensembl_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in ensembl.pep.all");
            //Console.WriteLine(gtf_seqs.Except(ensembl_seqs).Count().ToString() + " protein sequences only in gtf seqs");
            //Console.WriteLine();
            //Console.WriteLine("Comparison with Trembl sequences");
            //Console.WriteLine(trembl_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            //Console.WriteLine(trembl_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in Trembl");
            //Console.WriteLine(gtf_seqs.Except(trembl_seqs).Count().ToString() + " protein sequences only in gtf seqs");
            //Console.WriteLine();
            //Console.WriteLine("Comparison with SwissProt sequences");
            //Console.WriteLine(swiss_seqs.Intersect(gtf_seqs).Count().ToString() + " protein sequences in common");
            //Console.WriteLine(swiss_seqs.Except(gtf_seqs).Count().ToString() + " protein sequences only in Trembl");
            //Console.WriteLine(gtf_seqs.Except(swiss_seqs).Count().ToString() + " protein sequences only in gtf seqs");


            Console.ReadKey();

            ProteinDbWriter.WriteFastaDatabase(ensembl_proteins, @"C:\Users\antho\Documents\ProteoformDatabaseEngine\ensembl.fasta", " ");
            ProteinDbWriter.WriteFastaDatabase(uniq_pacbio_proteins, @"C:\Users\antho\Documents\ProteoformDatabaseEngine\pacbio.fasta", " ");

            //List<Protein> ensembl_seqs = ProteinDbLoader.LoadProteinFasta(@"C:\Users\antho\Documents\ProteoformDatabaseEngine\ensembl.fasta", true, DecoyType.None, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression).ToList();
            //List<Protein> pacbio_seqs = ProteinDbLoader.LoadProteinFasta(@"C:\Users\antho\Documents\ProteoformDatabaseEngine\pacbio.fasta", true, DecoyType.None, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression).ToList();

            Dictionary<string, List<Protein>> combined_seq_unique = new Dictionary<string, List<Protein>>();
            foreach (Protein p in ensembl_seqs.Concat(uniq_pacbio_proteins))
            {
                combined_seq_unique.TryGetValue(p.BaseSequence, out List<Protein> ps);
                if (ps == null) combined_seq_unique.Add(p.BaseSequence, new List<Protein> { p });
                else ps.Add(p);
            }
            List<Protein> proteins_combined = combined_seq_unique.Select(kv => kv.Value.First()).OrderBy(p => p.Accession).ToList();
            ProteinDbWriter.WriteFastaDatabase(proteins_combined, @"C:\Users\antho\Documents\ProteoformDatabaseEngine\ensembl_plus_pacbio.fasta", " ");
        }
    }
}
