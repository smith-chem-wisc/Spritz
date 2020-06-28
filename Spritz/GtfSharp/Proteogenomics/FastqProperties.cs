using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System;

namespace Proteogenomics
{
    public class FastqProperties
    {
        public FastqProperties(string fastqPath)
        {
            AssessReads(fastqPath);
        }

        public int ReadCount { get; set; }
        public double AverageReadLength { get; set; }

        private void AssessReads(string fastqPath)
        {
            int count = 0;
            List<int> readLengths = new List<int>();
            using (var stream = new FileStream(fastqPath, FileMode.Open, FileAccess.Read, FileShare.Read))
            {
                Stream fastaFileStream = fastqPath.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                StreamReader fastq = new StreamReader(fastaFileStream);

                while (true)
                {
                    string line = fastq.ReadLine();
                    if (line == null) { break; }
                    if (line.StartsWith("@"))
                    {
                        count++;
                        line = fastq.ReadLine();
                        if (line == null) { throw new NullReferenceException("FastqProperties error: strange file truncation."); }
                        readLengths.Add(line.Trim().Length);
                    }
                }
            }
            ReadCount = count;
            AverageReadLength = readLengths.Average();
        }
    }
}