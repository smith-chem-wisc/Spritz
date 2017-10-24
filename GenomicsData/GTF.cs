using System.Text.RegularExpressions;

namespace GenomicsData
{
    public class GTF
    {

        #region Private Fields

        private static Regex attribute_key = new Regex(@"([\w]+)"); // first instance of a word
        private static Regex attribute_value = new Regex(@"""([\w]+)"""); // anything inside the quotes

        #endregion Private Fields

        #region Public Method

        //public static GeneModel ReadGenomeFeatures(string gtf_location, Dictionary<string, Chromosome> chroms)
        //{

        //    GeneModel gene_model = new GeneModel(chroms, new List<Gene>());

        //    using (FileStream stream = new FileStream(gtf_location, FileMode.Open))
        //    {

        //        Stream vcfFileStream = gtf_location.EndsWith(".gz") ?
        //            (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
        //            stream;

        //        StreamReader gtf = new StreamReader(vcfFileStream);

        //        Gene gene = null;
        //        Transcript transcript = null;

        //        while (true)
        //        {
        //            string line = gtf.ReadLine();

        //            if (line == null)
        //                break;

        //            if (line.StartsWith("#"))
        //                continue;

        //            string[] fields = line.Split('\t');

        //            string seqname = fields[0];
        //            string source = fields[1];
        //            string feature = fields[2];
        //            int start = Convert.ToInt32(fields[3]);
        //            int end = Convert.ToInt32(fields[4]);
        //            string score = fields[5];
        //            string strand = fields[6];
        //            string frame = fields[7];
        //            Dictionary<string, string> attributes = new Dictionary<string, string>();
        //            foreach (string attrib in fields[8].Split(';').Where(x => x.Contains('"')))
        //            {
        //                string key = attribute_key.Match(attrib.TrimStart()).Groups[1].Value;
        //                string val = attribute_value.Match(attrib.TrimStart()).Groups[1].Value;
        //                if (!attributes.TryGetValue(key, out string x)) // sometimes there are two tags, so avoid adding twice
        //                    attributes.Add(key, val);
        //            }

        //            if (!chroms.TryGetValue(seqname, out Chromosome chrom))
        //                continue;

        //            bool has_gene_id = attributes.TryGetValue("gene_id", out string gene_id);
        //            bool has_gene_name = attributes.TryGetValue("gene_name", out string gene_name);
        //            bool has_gene_version = attributes.TryGetValue("gene_version", out string gene_version);
        //            bool has_gene_biotype = attributes.TryGetValue("gene_biotype", out string gene_biotype);
        //            bool has_transcript_id = attributes.TryGetValue("transcript_id", out string transcript_id);
        //            bool has_transcript_version = attributes.TryGetValue("transcript_version", out string transcript_version);
        //            bool has_transcript_biotype = attributes.TryGetValue("transcript_biotype", out string transcript_biotype);
        //            bool has_exon_id = attributes.TryGetValue("exon_id", out string exon_id);
        //            bool has_exon_version = attributes.TryGetValue("exon_version", out string exon_version);
        //            bool has_exon_number = attributes.TryGetValue("exon_number", out string exon_number);
        //            bool has_nearest_ref = attributes.TryGetValue("nearest_ref", out string nearest_ref); // Cufflinks
        //            bool has_class_code = attributes.TryGetValue("class_code", out string class_code); // Cufflinks
        //            int arbitrary_id = 1;

        //            switch (feature)
        //            {
        //                case "exon":
        //                    string exid = has_exon_id ? exon_id : has_exon_number && has_transcript_id ? transcript_id + "_" + exon_number : has_transcript_id ? transcript_id + "_" + arbitrary_id++.ToString() : arbitrary_id++.ToString();
        //                    if (gene == null || has_gene_id && gene_id != gene.ID)
        //                    {
        //                        gene = new Gene(gene_id, chrom, strand, start, end, gene_name, gene_biotype);
        //                        gene_model.genes.Add(gene);
        //                    }
        //                    if (transcript == null || has_transcript_id && transcript_id != transcript.ID)
        //                    {
        //                        transcript = new Transcript(transcript_id, chrom, strand, start, end, transcript_biotype, gene);
        //                        gene.transcripts.Add(transcript);
        //                    }
        //                    Exon exon = new Exon(exid, chrom, strand, start, end, transcript, gene);
        //                    gene.exons.Add(exon);
        //                    transcript.exons.Add(exon);
        //                    break;
        //                case "start_codon":
        //                    transcript.start_codon_start = start;
        //                    break;
        //                case "stop_codon":
        //                    transcript.stop_codon_start = start;
        //                    break;
        //                case "gene": // can be constructed from exons
        //                case "transcript": // can be constructed from exons
        //                case "CDS": // can be interpreted from exons & start/stop codons
        //                case "UTR": // can be interpreted from exons & start/stop codons
        //                    break;
        //            }
        //        }
        //    }
        //    return gene_model;
        //}

        #endregion Public Method

    }
}
