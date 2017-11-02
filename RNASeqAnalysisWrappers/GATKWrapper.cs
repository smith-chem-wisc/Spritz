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
        /// <summary>
        /// Groups and sorts reads, and marks duplicates using Picard Tools. Then, splits and trims reads with SplitNCigarReads.
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="bam"></param>
        /// <param name="new_bam"></param>
        public static void prepare_bam(string bin_directory, string bam, string genome_fasta, out string new_bam)
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
                groupsort_bam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + (sorted ? ".sorted.grouped.bam" : ".grouped.bam"));
                picard_command = "picard-tools AddOrReplaceReadGroups PU=platform  PL=illumina SM=sample LB=library I=" + WrapperUtility.convert_windows_path(bam) + " O=" + WrapperUtility.convert_windows_path(new_bam) + (sorted ? " SO=coordinate" : "");
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
                "java -jar GenomeAnalysisTK.jar -T SplitNCigarReads" + 
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) + 
                    " -I " + WrapperUtility.convert_windows_path(marked_duplicates_bam) + 
                    " -o " + WrapperUtility.convert_windows_path(split_trim_bam) + 
                    " -U ALLOW_N_CIGAR_READS",
                "rm " + WrapperUtility.convert_windows_path(marked_duplicates_bam)
            }).WaitForExit();
            new_bam = split_trim_bam;

            // run commands for marking duplicates and trimming reads
            //File.Delete(sorted_checkfile);
            //File.Delete(readgrouped_checkfile);
            //File.Delete(script_name);
            //File.Delete(script_name2);
            //File.Delete(script_name2);
        }

        /// <summary>
        /// Realigns indels for a given BAM file
        /// </summary>
        /// <param name="bin_directory"></param>
        /// <param name="genome_fasta"></param>
        /// <param name="bam"></param>
        /// <param name="known_sites_vcf"></param>
        /// <param name="new_bam"></param>
        public static void realign_indels(string bin_directory, string genome_fasta, string bam, string known_sites_vcf, out string new_bam)
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

                "java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator" +
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    " -known " + WrapperUtility.convert_windows_path(known_sites_vcf) +
                    " -o " +  WrapperUtility.convert_windows_path(realigner_table),

                "java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator" +
                    " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                    " -I " + WrapperUtility.convert_windows_path(bam) +
                    " -known " + WrapperUtility.convert_windows_path(known_sites_vcf) +
                    " -targetIntervals " +  WrapperUtility.convert_windows_path(realigner_table) +
                    " -o " + WrapperUtility.convert_windows_path(new_bam),
            }).WaitForExit();
        }

        // using BaseRecalibrator
        public static void base_recalibration(string bin_directory, string genome_fasta, string bam, string known_sites_vcf, out string recal_table_filepath)
        {
            recal_table_filepath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".recaltable");

            string arguments =
                " -T BaseRecalibrator" +
                " -R " + WrapperUtility.convert_windows_path(genome_fasta) +
                " -I " + WrapperUtility.convert_windows_path(bam) +
                " -knownSites " + WrapperUtility.convert_windows_path(known_sites_vcf) +
                " -o " + WrapperUtility.convert_windows_path(recal_table_filepath);

            string dictionary_path = Path.Combine(Path.GetDirectoryName(genome_fasta), Path.GetFileNameWithoutExtension(genome_fasta) + ".dict");
            string script_name = Path.Combine(bin_directory, "base_recalibration.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                !File.Exists(genome_fasta + ".fai") ? "samtools faidx " + WrapperUtility.convert_windows_path(genome_fasta) : "",
                !File.Exists(dictionary_path) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.convert_windows_path(genome_fasta) + " O=" + WrapperUtility.convert_windows_path(dictionary_path) : "",
                "java -jar GenomeAnalysisTK.jar" + arguments
            }).WaitForExit();
            //File.Delete(script_name);
        }

        // using HaplotypeCaller on each BAM file individually
        public static Process variant_calling()
        {
            return null;
        }

        // using GenotypeGVCFs on all VCF files together
        public static Process genotype()
        {
            return null;
        }

        // using VariantRecalibrator, ApplyRecalibration, and then VariantFiltration
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
            if (File.Exists(Path.Combine(target_directory, known_sites_filename))) return;
            string script_path = Path.Combine(bin_directory, "download_known_variants.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(target_directory),
                "wget " + target_file,
                "tar -xvf " + target_file,
                "rm " + target_file
            }).WaitForExit();
            File.Delete(script_path);
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
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
