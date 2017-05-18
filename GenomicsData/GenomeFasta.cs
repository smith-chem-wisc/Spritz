using Genomics;
using Ionic.Zlib;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace GenomicsData
{
    public class GenomeFasta
    {

        #region Public Methods

        public static Dictionary<string, Chromosome> ReadGenomeFasta(string genomeFastaLocation)
        {
            string chrom_name = null;
            int unique_identifier = 1;
            HashSet<string> unique_names = new HashSet<string>();
            Dictionary<string, Chromosome> chromosomes = new Dictionary<string, Chromosome>();
            Regex substituteWhitespace = new Regex(@"\s+");

            using (FileStream stream = new FileStream(genomeFastaLocation, FileMode.Open))
            {
                Stream fastaFileStream = genomeFastaLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                StringBuilder sb = null;
                StreamReader fasta = new StreamReader(fastaFileStream);

                while (true)
                {
                    string line = fasta.ReadLine();

                    if (line.StartsWith(">"))
                    {
                        chrom_name = substituteWhitespace.Split(line.Substring(1).TrimEnd())[0];
                        sb = new StringBuilder();
                    }

                    else if (sb != null)
                    {
                        sb.Append(line.Trim());
                    }

                    if ((fasta.Peek() == '>' || fasta.Peek() == -1) & chrom_name != null && sb != null)
                    {
                        string sequence = substituteWhitespace.Replace(sb.ToString(), "");
                        while (unique_names.Contains(chrom_name))
                        {
                            chrom_name += "_" + unique_identifier.ToString();
                            unique_identifier++;
                        }
                        unique_names.Add(chrom_name);
                        Chromosome chrom = new Chromosome(chrom_name, sequence);
                        chromosomes.Add(chrom_name, chrom);

                        chrom_name = null;
                    }

                    if (fasta.Peek() == -1)
                    {
                        break;
                    }
                }
            }

            return chromosomes;
        }

        #endregion Public Methods
    }
}
