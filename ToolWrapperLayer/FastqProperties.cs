using System.IO;
using System.IO.Compression;

namespace ToolWrapperLayer
{
    public class FastqProperties
    {
        public FastqProperties(string fastqPath)
        {
            ReadCount = CountReads(fastqPath);
        }

        public int ReadCount { get; set; }

        private int CountReads(string fastqPath)
        {
            int count = 0;
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
                    if (line.StartsWith("@")) { count++; }
                }
            }
            return count;
        }
    }
}