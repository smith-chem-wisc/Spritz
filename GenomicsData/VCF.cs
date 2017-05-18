using Genomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using Ionic.Zlib;

namespace GenomicsData
{
    /// <summary>
    /// Variant Call Format, https://samtools.github.io/hts-specs/VCFv4.2.pdf
    /// </summary>
    public class VCF
    {

        #region Private Fields

        Regex gt_allele = new Regex(@"(\d+)");
        Regex gt_phase = new Regex(@"([|/])");

        #endregion Private Fields

        #region Public Methods

        public static List<Sample> ReadVCF(string vcfLocation, Dictionary<string, Chromosome> chroms)
        {
            List<Sample> samples = new List<Sample>();
            using (FileStream stream = new FileStream(vcfLocation, FileMode.Open))
            {

                Stream vcfFileStream = vcfLocation.EndsWith(".gz") ?
                    (Stream)(new GZipStream(stream, CompressionMode.Decompress)) :
                    stream;

                StreamReader vcf = new StreamReader(vcfFileStream);

                while (true)
                {
                    string line = vcf.ReadLine();

                    if (line == null)
                        break;

                    string vcf_format = "";
                    if (line.StartsWith("##")) //metainfo
                    {
                        if (line.StartsWith("##fileformat"))
                            vcf_format = line.Split('=')[1];
                        continue;
                    }

                    string[] fields = line.Split('\t');
                    if (line.StartsWith("#")) //header
                    {
                        samples = Enumerable.Range(9, fields.Length - 9).Select(i => new Sample(fields[i])).ToList();
                        continue;
                    }

                    if (samples.Count <= 0) // no header or no sample information
                        samples.Add(new Sample(""));

                    string chrom_name = fields[0];
                    int oneBasedPosition = Convert.ToInt32(fields[1]);
                    string id = fields[2];
                    string reference = fields[3];
                    string[] alternate = fields[4].Split(',');
                    double qual = Convert.ToDouble(fields[5]);
                    string filter = fields[6];
                    Dictionary<string, string> info = fields[7].Split(';').ToDictionary(x => x.Split('=').First(), x => x.Split('=').Last());
                    string[] format = fields[8].Split(':');
                    Dictionary<string, string> sample_specs = Enumerable.Range(9, fields.Length - 9).ToDictionary(i => samples[i - 9].name, i => fields[i]);

                    if (!chroms.TryGetValue(chrom_name, out Chromosome chrom))
                        continue;

                    //bool last_phase = false;
                    //LocalHaplotype local_haplotype = null;
                    foreach (Sample s in samples)
                    {
                        bool sample_specified = sample_specs.TryGetValue(s.name, out string sample_spec);
                        int allele_num = 1;

                        //int allele_num = sample_specified && sample_spec.Contains("GT") ?
                        //  Convert.ToInt32(gt_allele.Match(sample_spec).Groups[1].Value) :
                        //  1;

                        //bool in_phase = sample_specified && sample_spec.Contains("HP") && gt_phase.Match(sample_spec).Groups[1].Value == "|"; //need to fix this for the GATK output
                        //if (!last_phase && in_phase)
                        //{
                        //    local_haplotype = new LocalHaplotype(chrom, pos);
                        //    s.local_haplotypes.Add(local_haplotype);
                        //}

                        SequenceVariant seqvar;

                        if (reference.Length > alternate.Length)
                            seqvar = new Deletion(chrom, oneBasedPosition, id, reference, alternate[allele_num - 1], qual, filter, info);

                        else if (reference.Length < alternate.Length)
                            seqvar = new Insertion(chrom, oneBasedPosition, id, reference, alternate[allele_num - 1], qual, filter, info);

                        else
                            seqvar = new SNV(chrom, oneBasedPosition, id, reference, alternate[allele_num - 1], qual, filter, info);

                        s.sequence_variants.Add(seqvar);
                        //if (local_haplotype != null && in_phase) local_haplotype.add()
                    }
                }
            }
            return samples;
        }

        #endregion Public Methods

    }
}
