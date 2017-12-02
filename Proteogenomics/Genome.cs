using Bio;
using Bio.IO.FastA;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System;

namespace Proteogenomics
{
    public class Genome
    {

        #region Public Properties

        public List<ISequence> Chromosomes { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public Genome(string genomeFastaLocation)
        {
            Chromosomes = new FastAParser().Parse(genomeFastaLocation).ToList();
        }

        #endregion Public Constructor


        #region Public Method

        public List<ISequence> KaryotypicOrder()
        {
            ISequence[] orderedChromosomes = new ISequence[Chromosomes.Count];
            bool ucsc = Chromosomes[0].ID.StartsWith("c");
            int i = 0;
            foreach (int chr in Enumerable.Range(1, 22))
            {
                ISequence s = Chromosomes.FirstOrDefault(x => x.ID.Split(' ')[0] == (ucsc ? "chr" + chr : chr.ToString()));
                if (s != null) orderedChromosomes[i++] = s;
            }
            ISequence seqx = Chromosomes.FirstOrDefault(x => x.ID.Split(' ')[0] == (ucsc ? "chrX" : "X"));
            if (seqx != null) orderedChromosomes[i++] = seqx;
            ISequence seqy = Chromosomes.FirstOrDefault(x => x.ID.Split(' ')[0] == (ucsc ? "chrY" : "Y"));
            if (seqy != null) orderedChromosomes[i++] = seqy;
            ISequence seqm = Chromosomes.FirstOrDefault(x => x.ID.Split(' ')[0] == (ucsc ? "chrM" : "MT"));
            if (seqm != null) orderedChromosomes[i++] = seqm;

            List<ISequence> gl = Chromosomes.Where(x => x.ID.Split(' ').Contains("GL")).ToList();
            foreach (var g in gl)
            {
                orderedChromosomes[i++] = g;
            }

            List<ISequence> ki = Chromosomes.Where(x => x.ID.Split(' ').Contains("KI")).ToList();
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

        public static void WriteFasta(IEnumerable<ISequence> sequences, string filePath)
        {
            FastAFormatter formatter = new FastAFormatter();
            using (FileStream stream = File.Create(filePath))
                formatter.Format(stream, sequences);
            using (StreamReader reader = new StreamReader(filePath))
            using (StreamWriter writer = new StreamWriter(filePath + ".tmp"))
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null) break;
                    writer.Write(line + '\n');
                }
            File.Delete(filePath);
            File.Move(filePath + ".tmp", filePath);
        }

        #endregion Public Method
    }
}
