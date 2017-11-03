using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace RNASeqAnalysisWrappers
{
    public class GATKWrapper
    {
        public static string gatk = "java -Xmx20G -jar GenomeAnalysisTK.jar";

        public static void subset_bam(string bin_directory, int threads, string bam, string genome_fasta, string genome_region, string output_bam)
        {
            string script_name = Path.Combine(bin_directory, "subset_bam.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                gatk +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    " -o " + WrapperUtility.convert_windows_path(output_bam) +
                    " -L " + genome_region,
            }).WaitForExit();
        }

        /// <summary>
        /// Groups and sorts reads, and marks duplicates using Picard Tools. Then, splits and trims reads with SplitNCigarReads.
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="bam"></param>
        /// <param name="new_bam"></param>
        public static void prepare_bam(string bin_directory, int threads, string bam, string genome_fasta, out string new_bam)
        {
            new_bam = bam;

            // check if sorted and grouped
            string sorted_checkfile = "header_sorted.txt";
            string readgrouped_checkfile = "header_readgrouped.txt";
            string script_name = Path.Combine(bin_directory, "check_sorted.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "samtools view -H " + WrapperUtility.convert_windows_path(bam) + " | grep SO:coordinate > " + sorted_checkfile,
                "samtools view -H " + WrapperUtility.convert_windows_path(bam) + " | grep '^@RG' > " + readgrouped_checkfile,
            }).WaitForExit();

            bool sorted = new FileInfo(Path.Combine(bin_directory, sorted_checkfile)).Length > 0;
            bool grouped = new FileInfo(Path.Combine(bin_directory, readgrouped_checkfile)).Length > 0;

            // run commands for grouping or sorting
            string picard_command;
            string groupsort_bam = new_bam;
            if (grouped && sorted) return;
            else if (!grouped)
            {
                groupsort_bam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + (!sorted ? ".sorted.grouped.bam" : ".grouped.bam"));
                picard_command = "picard-tools AddOrReplaceReadGroups PU=platform  PL=illumina SM=sample LB=library I=" + WrapperUtility.convert_windows_path(bam) + " O=" + WrapperUtility.convert_windows_path(groupsort_bam) + (!sorted ? " SO=coordinate" : "");
            }
            else // grouped, but not sorted
            {
                groupsort_bam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".sorted.bam") ;
                picard_command = "picard-tools SortSam SO=coordinate I=" + WrapperUtility.convert_windows_path(bam) + " O=" + WrapperUtility.convert_windows_path(groupsort_bam);
            }

            string marked_duplicates_bam = Path.Combine(Path.GetDirectoryName(groupsort_bam), Path.GetFileNameWithoutExtension(groupsort_bam) + ".marked.bam");
            string marked_duplicate_metrics = Path.Combine(Path.GetDirectoryName(groupsort_bam), Path.GetFileNameWithoutExtension(groupsort_bam) + ".marked.metrics");
            string dictionary_path = Path.Combine(Path.GetDirectoryName(genome_fasta), Path.GetFileNameWithoutExtension(genome_fasta) + ".dict");
            string split_trim_bam = Path.Combine(Path.GetDirectoryName(marked_duplicates_bam), Path.GetFileNameWithoutExtension(marked_duplicates_bam) + ".split.bam");
            string script_name2 = Path.Combine(bin_directory, "picard.bash");
            WrapperUtility.generate_and_run_script(script_name2, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                picard_command,
                "picard-tools MarkDuplicates I=" + WrapperUtility.convert_windows_path(groupsort_bam) +
                    " O=" + WrapperUtility.convert_windows_path(marked_duplicates_bam) +
                    " M=" + WrapperUtility.convert_windows_path(marked_duplicate_metrics) +
                    " AS=true", // assume sorted
                "rm " + WrapperUtility.convert_windows_path(groupsort_bam), // conserve space

                !File.Exists(genome_fasta + ".fai") ? "samtools faidx " + WrapperUtility.convert_windows_path(genome_fasta) : "",
                !File.Exists(dictionary_path) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.convert_windows_path(genome_fasta) + " O=" + WrapperUtility.convert_windows_path(dictionary_path) : "",

                "samtools index " + WrapperUtility.convert_windows_path(marked_duplicates_bam),
                gatk +
                    " -T SplitNCigarReads" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) + 
                    " -I " + WrapperUtility.convert_windows_path(marked_duplicates_bam) + 
                    " -o " + WrapperUtility.convert_windows_path(split_trim_bam) + 
                    " -U ALLOW_N_CIGAR_READS -fixMisencodedQuals", // STAR apparently misencodes quality scores; this just subtracts 31 from all quality scores... might not need that flag for tophat
                "rm " + WrapperUtility.convert_windows_path(marked_duplicates_bam),
                "rm " + WrapperUtility.convert_windows_path(Path.Combine(Path.GetDirectoryName(marked_duplicates_bam), Path.GetFileNameWithoutExtension(marked_duplicates_bam) + ".bai")),
            }).WaitForExit();
            new_bam = split_trim_bam;

            // run commands for marking duplicates and trimming reads
            File.Delete(sorted_checkfile);
            File.Delete(readgrouped_checkfile);
            File.Delete(script_name);
            File.Delete(script_name2);
            File.Delete(script_name2);
        }

        /// <summary>
        /// Realigns indels for a given BAM file
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="genome_fasta"></param>
        /// <param name="bam"></param>
        /// <param name="known_sites_vcf"></param>
        /// <param name="new_bam"></param>
        public static void realign_indels(string bin_directory, int threads, string genome_fasta, string bam, out string new_bam, string known_sites_vcf = "")
        {
            string dictionary_path = Path.Combine(Path.GetDirectoryName(genome_fasta), Path.GetFileNameWithoutExtension(genome_fasta) + ".dict");
            string realigner_table = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".forIndelRealigner.intervals");
            new_bam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".realigned.bam");
            string script_name = Path.Combine(bin_directory, "realign_indels.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                !File.Exists(genome_fasta + ".fai") ? "samtools faidx " + WrapperUtility.convert_windows_path(genome_fasta) : "",
                !File.Exists(dictionary_path) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.convert_windows_path(genome_fasta) + " O=" + WrapperUtility.convert_windows_path(dictionary_path) : "",

                gatk + 
                    " -T RealignerTargetCreator" +
                    " --num_threads " + threads.ToString() +  
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    (known_sites_vcf != "" ? " -known " + WrapperUtility.convert_windows_path(known_sites_vcf) : "") +
                    " -o " +  WrapperUtility.convert_windows_path(realigner_table),

                gatk + 
                    " -T IndelRealigner" +
                    //" --num_threads " + threads.ToString() + // this tool can't do threaded analysis
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    (known_sites_vcf != "" ? " -known " + WrapperUtility.convert_windows_path(known_sites_vcf) : "") +
                    " -targetIntervals " +  WrapperUtility.convert_windows_path(realigner_table) +
                    " -o " + WrapperUtility.convert_windows_path(new_bam),
            }).WaitForExit();
        }

        /// <summary>
        /// Creates recalibration table for base calls. 
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="genome_fasta"></param>
        /// <param name="bam"></param>
        /// <param name="recal_table_filepath"></param>
        /// <param name="known_sites_vcf"></param>
        public static void base_recalibration(string bin_directory, string genome_fasta, string bam, out string recal_table_filepath, string known_sites_vcf) 
        {
            recal_table_filepath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".recaltable");

            string dictionary_path = Path.Combine(Path.GetDirectoryName(genome_fasta), Path.GetFileNameWithoutExtension(genome_fasta) + ".dict");
            string script_name = Path.Combine(bin_directory, "base_recalibration.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                !File.Exists(genome_fasta + ".fai") ? "samtools faidx " + WrapperUtility.convert_windows_path(genome_fasta) : "",
                !File.Exists(dictionary_path) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.convert_windows_path(genome_fasta) + " O=" + WrapperUtility.convert_windows_path(dictionary_path) : "",

                gatk + 
                    " -T BaseRecalibrator" +
                    //" --num_threads " + threads.ToString() + // doesn't support threaded runs
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    (known_sites_vcf != "" ? " -knownSites " + WrapperUtility.convert_windows_path(known_sites_vcf) : "") +
                    " -o " + WrapperUtility.convert_windows_path(recal_table_filepath),
            }).WaitForExit();
            File.Delete(script_name);
        }

        /// <summary>
        /// HaplotypeCaller for calling variants on each RNA-Seq BAM file individually
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="threads"></param>
        /// <param name="genome_fasta"></param>
        /// <param name="bam"></param>
        /// <param name=""></param>
        public static void variant_calling(string bin_directory, int threads, string genome_fasta, string bam, string dbsnp, out string new_vcf)
        {
            new_vcf = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".vcf");
            string dictionary_path = Path.Combine(Path.GetDirectoryName(genome_fasta), Path.GetFileNameWithoutExtension(genome_fasta) + ".dict");
            string script_name = Path.Combine(bin_directory, "variant_calling.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                !File.Exists(genome_fasta + ".fai") ? "samtools faidx " + WrapperUtility.convert_windows_path(genome_fasta) : "",
                !File.Exists(dictionary_path) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.convert_windows_path(genome_fasta) + " O=" + WrapperUtility.convert_windows_path(dictionary_path) : "",

                gatk +
                    " -T HaplotypeCaller" +
                    " -nct " + threads.ToString() +
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    " --standard_min_confidence_threshold_for_calling 20" +
                    " --dbsnp " + WrapperUtility.convert_windows_path(dbsnp) +
                    " -o " + WrapperUtility.convert_windows_path(new_vcf)
            }).WaitForExit();
        }

        // using GenotypeGVCFs on all VCF files together
        public static Process genotype()
        {
            return null;
        }

        // using VariantRecalibrator, ApplyRecalibration, and then VariantFiltration -- apparently mostly for DNAseq
        public static Process filter_variants()
        {
            return null;
        }

        private static string all_grch37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz";
        private static string common_grch37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/common_all_20170710.vcf.gz";
        private static string all_grch38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz";
        private static string common_grch38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz";

        public static void download_known_sites(string bin_directory, string target_directory, bool common_only, bool GRCh37, bool GRCh38, out string known_sites_filename)
        {
            known_sites_filename = "";

            if (!GRCh37 && !GRCh38) return;

            string target_file = common_only ?
                GRCh37 ? common_grch37 : common_grch38 :
                GRCh37 ? all_grch37 : all_grch38;

            known_sites_filename = target_file.Split('/').Last();
            known_sites_filename = known_sites_filename.Substring(0, known_sites_filename.Length - 3);
            string new_known_sites = Path.GetFileNameWithoutExtension(known_sites_filename) + ".ensembl.vcf";

            if (File.Exists(Path.Combine(target_directory, new_known_sites)))
            {
                known_sites_filename = new_known_sites;
                return;
            }

            string script_path = Path.Combine(bin_directory, "download_known_variants.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(target_directory),
                "wget " + target_file,
                "tar -xvf " + target_file,
                "rm " + target_file
            }).WaitForExit();
            File.Delete(script_path);

            Dictionary<string, string> chrom_mappings = File.ReadAllLines(GRCh37 ? Path.Combine(bin_directory, "ChromosomeMappings", "GRCh37_UCSC2ensembl.txt") : Path.Combine(bin_directory, "ChromosomeMappings", "GRCh38_UCSC2ensembl.txt"))
                .ToDictionary(line => line.Split('\t')[0], line => line.Split('\t')[1]);

            using (StreamReader reader = new StreamReader(Path.Combine(target_directory, known_sites_filename)))
            using (StreamWriter writer = new StreamWriter(Path.Combine(target_directory, new_known_sites)))
            {
                while (true)
                {
                    string a = reader.ReadLine();
                    if (a == null) break;
                    string[] line = a.Split('\t');
                    if (chrom_mappings.TryGetValue(line[0], out string chr)) line[0] = chr;
                    writer.Write(String.Join("\t", line) + '\n');
                }
            }
            known_sites_filename = new_known_sites;
            //File.Delete(Path.Combine(target_directory, known_sites_filename));
        }

        public static void install(string current_directory)
        {
            if (Directory.Exists(Path.Combine(current_directory, "GenomeAnalysisTK.jar"))) return;
            string script_path = Path.Combine(current_directory, "install_gatk.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                @"wget https://software.broadinstitute.org/gatk/download/auth?package=GATK",
                "tar -jxvf GenomeAnalysisTK-3.8-0.tar.bz2",
                "rm GenomeAnalysisTK-3.8-0.tar.bz2",
                @"mv GenomeAnalysisTK-*/GenomeAnalysisTK.jar .",
                @"rm -r GenomeAnalysisTK-*",
                "git clone https://github.com/dpryan79/ChromosomeMappings.git",
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
