using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Genomics;

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

        #region Public Properties

        public string format { get; set; }
        public List<Sample> samples { get; set; }

        #endregion Public Properties

        #region Public Constructor

        public VCF(string filepath)
        {
            using (StreamReader reader = new StreamReader(filepath))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line.StartsWith("##")) //metainfo
                    {
                        if (line.StartsWith("##fileformat"))
                            this.format = line.Split('=')[1];
                        continue;
                    }

                    string[] fields = line.Split('\t');
                    if (line.StartsWith("#")) //header
                    {
                        samples = Enumerable.Range(9, fields.Length - 9).Select(i => new Sample(fields[i])).ToList();
                        if (samples.Count <= 0)
                            samples.Add(new Sample(""));
                        continue;
                    }

                    string chrom_name = fields[0];
                    int pos = Convert.ToInt32(fields[1]);
                    string id = fields[2];
                    string reference = fields[3];
                    string[] alternate = fields[4].Split(',');
                    double qual = Convert.ToDouble(fields[5]);
                    string filter = fields[6];
                    Dictionary<string, string> info = fields[7].Split(';').ToDictionary(x => x.Split('=').First(), x => x.Split('=').Last());
                    string[] format = fields[8].Split(':');
                    Dictionary<string, string> sample_specs = Enumerable.Range(9, fields.Length - 9).ToDictionary(i => samples[i].name, i => fields[i]);

                    //bool last_phase = false;
                    //LocalHaplotype local_haplotype = null;
                    foreach (Sample s in samples)
                    {
                        if (!s.chroms.TryGetValue(chrom_name, out Chromosome chrom))
                        {
                            chrom = new Chromosome(chrom_name);
                            s.chroms.Add(chrom_name, chrom);
                        }

                        bool sample_specified = sample_specs.TryGetValue(s.name, out string sample_spec);
                        int allele_num = 1;

                        //int allele_num = sample_specified && sample_spec.Contains("GT") ?
                        //  Convert.ToInt32(gt_allele.Match(sample_spec).Value) :
                        //  1;

                        //bool in_phase = sample_specified && sample_spec.Contains("HP") && gt_phase.Match(sample_spec).Value == "|"; //need to fix this for the GATK output
                        //if (!last_phase && in_phase)
                        //{
                        //    local_haplotype = new LocalHaplotype(chrom, pos);
                        //    s.local_haplotypes.Add(local_haplotype);
                        //}

                        SequenceVariant seqvar;
                        if (reference.Length > alternate.Length)
                        {
                            seqvar = new Deletion(chrom, pos, id, reference, alternate[allele_num - 1], qual, filter, info);
                        }
                        else if (reference.Length < alternate.Length)
                        {
                            seqvar = new Insertion(chrom, pos, id, reference, alternate[allele_num - 1], qual, filter, info);
                        }
                        else
                        {
                            seqvar = new SNV(chrom, pos, id, reference, alternate[allele_num - 1], qual, filter, info);
                        }

                        s.sequence_variants.Add(seqvar);
                        //if (local_haplotype != null && in_phase) local_haplotype.add()
                    }
                }
            }
        }

        #endregion Public Constructor

    }
}
