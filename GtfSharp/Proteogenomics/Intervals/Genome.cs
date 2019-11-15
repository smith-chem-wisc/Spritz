using Bio;
using Bio.IO.FastA;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Proteogenomics
{
    public class Genome
        : Interval
    {
        /// <summary>
        /// Parses the chromosomes listed in a genome fasta.
        /// </summary>
        /// <param name="genomeFastaLocation"></param>
        public Genome(string genomeFastaLocation)
        {
            Chromosomes = new FastAParser().Parse(genomeFastaLocation).Select(s => new Chromosome(s, this)).ToList();
        }

        /// <summary>
        /// Chromosomes contained in this genome.
        /// </summary>
        public List<Chromosome> Chromosomes { get; set; }

        /// <summary>
        /// Orders the chromosomes in this genome by karyotypic order, i.e. by chromosome size.
        /// </summary>
        /// <param name="fastaHeaderDelimeter"></param>
        /// <returns></returns>
        public List<Chromosome> KaryotypicOrder()
        {
            Chromosome[] orderedChromosomes = new Chromosome[Chromosomes.Count];
            bool ucsc = Chromosomes[0].FriendlyName.StartsWith("c");
            int i = 0;
            foreach (int chr in Enumerable.Range(1, 22))
            {
                Chromosome s = Chromosomes.FirstOrDefault(x => x.FriendlyName == (ucsc ? "chr" + chr : chr.ToString()));
                if (s != null) { orderedChromosomes[i++] = s; }
            }
            Chromosome seqx = Chromosomes.FirstOrDefault(x => x.FriendlyName == (ucsc ? "chrX" : "X"));
            if (seqx != null) { orderedChromosomes[i++] = seqx; }
            Chromosome seqy = Chromosomes.FirstOrDefault(x => x.FriendlyName == (ucsc ? "chrY" : "Y"));
            if (seqy != null) { orderedChromosomes[i++] = seqy; }
            Chromosome seqm = Chromosomes.FirstOrDefault(x => x.FriendlyName == (ucsc ? "chrM" : "MT"));
            if (seqm != null) { orderedChromosomes[i++] = seqm; }

            List<Chromosome> gl = Chromosomes.Where(x => x.FriendlyName.Contains("GL")).ToList();
            foreach (var g in gl)
            {
                orderedChromosomes[i++] = g;
            }

            List<Chromosome> ki = Chromosomes.Where(x => x.FriendlyName.Contains("KI")).ToList();
            foreach (var k in ki)
            {
                orderedChromosomes[i++] = k;
            }

            foreach (var x in Chromosomes.Except(orderedChromosomes))
            {
                orderedChromosomes[i++] = x;
            }
            return orderedChromosomes.ToList();
        }

        /// <summary>
        /// Checks whether the chromosomes are listed in karyotypic order, i.e. by chromosome size.
        /// </summary>
        /// <param name="fastaHeaderDelimeter"></param>
        /// <returns></returns>
        public bool IsKaryotypic()
        {
            bool ucsc = Chromosomes[0].FriendlyName.StartsWith("c");
            int i = 0;
            List<string> ids = Chromosomes.Select(x => x.FriendlyName).ToList();
            List<string> names = new List<string>();
            foreach (string chr in Enumerable.Range(1, 22).Select(x => x.ToString()).Concat(new string[] { "X", "Y", "M" }))
            {
                string name = ucsc ? "chr" + chr : chr.ToString() + (chr == "M" ? "T" : "");
                names.Add(name);
                int s = ids.IndexOf(name);
                if (s > 0)
                {
                    i = s;
                }
                if (s > 0 && s <= i)
                {
                    return false;
                }
            }
            foreach (string chr in ids.Except(names))
            {
                string name = ucsc ? "chr" + chr : chr.ToString() + (chr == "M" ? "T" : "");
                int s = ids.IndexOf(name);
                if (s > 0 && s <= i)
                {
                    return false;
                }
            }
            return true;
        }

        public static void WriteFasta(IEnumerable<ISequence> sequences, string filePath)
        {
            FastAFormatter formatter = new FastAFormatter();
            using (FileStream stream = File.Create(filePath))
                formatter.Format(stream, sequences);
            using (StreamReader reader = new StreamReader(filePath))
            using (StreamWriter writer = new StreamWriter(filePath + ".tmp"))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null) { break; }
                    writer.Write(line + '\n');
                }
            }
            File.Delete(filePath);
            File.Move(filePath + ".tmp", filePath);
        }
    }
}