using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace RNASeqAnalysisWrappers
{
    public class GATKWrapper
    {
        public static string GATK = "java -Xmx20G -jar GenomeAnalysisTK.jar";

        public static void SubsetBam(string binDirectory, int threads, string bam, string genomeFasta, string genomeRegion, string outputBam)
        {
            string script_name = Path.Combine(binDirectory, "scripts", "subset_bam.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GATK +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(outputBam) +
                    " -L " + genomeRegion,
            }).WaitForExit();
        }

        /// <summary>
        /// Groups and sorts reads, and marks duplicates using Picard Tools. Then, splits and trims reads with SplitNCigarReads.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="bam"></param>
        /// <param name="newBam"></param>
        public static void PrepareBam(string binDirectory, int threads, string bam, string genomeFasta, out string newBam)
        {
            newBam = bam;

            // check if sorted and grouped
            string sortedCheckfile = "header_sorted.txt";
            string readGroupedCheckfile = "header_readgrouped.txt";
            string scriptPath = Path.Combine(binDirectory, "scripts", "check_sorted.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep SO:coordinate > " + sortedCheckfile,
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep '^@RG' > " + readGroupedCheckfile,
            }).WaitForExit();

            bool sorted = new FileInfo(Path.Combine(binDirectory, sortedCheckfile)).Length > 0;
            bool grouped = new FileInfo(Path.Combine(binDirectory, readGroupedCheckfile)).Length > 0;

            // run commands for grouping or sorting
            string picardCommand;
            string groupSortBam = newBam;
            if (grouped && sorted) return;
            else if (!grouped)
            {
                groupSortBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + (!sorted ? ".sorted.grouped.bam" : ".grouped.bam"));
                picardCommand = "picard-tools AddOrReplaceReadGroups PU=platform  PL=illumina SM=sample LB=library I=" + WrapperUtility.ConvertWindowsPath(bam) + " O=" + WrapperUtility.ConvertWindowsPath(groupSortBam) + (!sorted ? " SO=coordinate" : "");
            }
            else // grouped, but not sorted
            {
                groupSortBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".sorted.bam") ;
                picardCommand = "picard-tools SortSam SO=coordinate I=" + WrapperUtility.ConvertWindowsPath(bam) + " O=" + WrapperUtility.ConvertWindowsPath(groupSortBam);
            }

            string markedDuplicatesBam = Path.Combine(Path.GetDirectoryName(groupSortBam), Path.GetFileNameWithoutExtension(groupSortBam) + ".marked.bam");
            string markedDuplicateMetrics = Path.Combine(Path.GetDirectoryName(groupSortBam), Path.GetFileNameWithoutExtension(groupSortBam) + ".marked.metrics");
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".dict");
            string splitTrimBam = Path.Combine(Path.GetDirectoryName(markedDuplicatesBam), Path.GetFileNameWithoutExtension(markedDuplicatesBam) + ".split.bam");
            string mapQReassigned = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".mapqfixed.bam");
            string scriptName2 = Path.Combine(binDirectory, "scripts", "picard.bash");
            WrapperUtility.GenerateAndRunScript(scriptName2, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                picardCommand,
                "picard-tools MarkDuplicates I=" + WrapperUtility.ConvertWindowsPath(groupSortBam) +
                    " O=" + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) +
                    " M=" + WrapperUtility.ConvertWindowsPath(markedDuplicateMetrics) +
                    " AS=true", // assume sorted
                "rm " + WrapperUtility.ConvertWindowsPath(groupSortBam), // conserve space

                !File.Exists(genomeFasta + ".fai") ? "samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFasta) : "",
                !File.Exists(dictionaryPath) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.ConvertWindowsPath(genomeFasta) + " O=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) : "",

                "samtools index " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam),
                GATK +
                    " -T SplitNCigarReads" +
                    //" --num_threads " + threads.ToString() + // not supported
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " -U ALLOW_N_CIGAR_READS -fixMisencodedQuals", // STAR apparently misencodes quality scores; this just subtracts 31 from all quality scores... might not need that flag for tophat
                "rm " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam),
                "rm " + WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(markedDuplicatesBam), Path.GetFileNameWithoutExtension(markedDuplicatesBam) + ".bai")),
                GATK +
                    " -T PrintReads" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(mapQReassigned) +
                    " -rf ReassignMappingQuality", // default mapping quality is 60; required for RNA-Seq aligners
                "samtools index " + WrapperUtility.ConvertWindowsPath(mapQReassigned),
                "rm " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam),
            }).WaitForExit();
            newBam = mapQReassigned;

            // run commands for marking duplicates and trimming reads
            File.Delete(sortedCheckfile);
            File.Delete(readGroupedCheckfile);
            //File.Delete(scriptName);
            //File.Delete(scriptName2);
        }

        /// <summary>
        /// Realigns indels for a given BAM file
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="bam"></param>
        /// <param name="knownSitesVcf"></param>
        /// <param name="newBam"></param>
        public static void RealignIndels(string binDirectory, int threads, string genomeFasta, string bam, out string newBam, string knownSitesVcf = "")
        {
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".dict");
            string realignerTable = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".forIndelRealigner.intervals");
            newBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".realigned.bam");
            string script_name = Path.Combine(binDirectory, "scripts", "realign_indels.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                !File.Exists(genomeFasta + ".fai") ? "samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFasta) : "",
                !File.Exists(dictionaryPath) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.ConvertWindowsPath(genomeFasta) + " O=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) : "",

                GATK + 
                    " -T RealignerTargetCreator" +
                    " --num_threads " + threads.ToString() +  
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -o " +  WrapperUtility.ConvertWindowsPath(realignerTable),

                GATK + 
                    " -T IndelRealigner" +
                    //" --num_threads " + threads.ToString() + // this tool can't do threaded analysis
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -targetIntervals " +  WrapperUtility.ConvertWindowsPath(realignerTable) +
                    " -o " + WrapperUtility.ConvertWindowsPath(newBam),
            }).WaitForExit();
        }

        /// <summary>
        /// Creates recalibration table for base calls. 
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="bam"></param>
        /// <param name="recalibrationTablePath"></param>
        /// <param name="knownSitesVcf"></param>
        public static void BaseRecalibration(string binDirectory, string genomeFasta, string bam, out string recalibrationTablePath, string knownSitesVcf) 
        {
            recalibrationTablePath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".recaltable");

            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".dict");
            string scriptPath = Path.Combine(binDirectory, "scripts", "base_recalibration.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                !File.Exists(genomeFasta + ".fai") ? "samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFasta) : "",
                !File.Exists(dictionaryPath) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.ConvertWindowsPath(genomeFasta) + " O=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) : "",

                GATK + 
                    " -T BaseRecalibrator" +
                    //" --num_threads " + threads.ToString() + // doesn't support threaded runs
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -knownSites " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -o " + WrapperUtility.ConvertWindowsPath(recalibrationTablePath),
            }).WaitForExit();
        }

        /// <summary>
        /// HaplotypeCaller for calling variants on each RNA-Seq BAM file individually
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="bam"></param>
        /// <param name=""></param>
        public static void VariantCalling(string binDirectory, int threads, string genomeFasta, string bam, string dbsnp, out string newVcf)
        {
            newVcf = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".vcf");

            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".dict");
            string scriptName = Path.Combine(binDirectory, "scripts", "variant_calling.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                !File.Exists(genomeFasta + ".fai") ? "samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFasta) : "",
                !File.Exists(dictionaryPath) ? "picard-tools CreateSequenceDictionary R=" + WrapperUtility.ConvertWindowsPath(genomeFasta) + " O=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) : "",

                GATK +
                    " -T HaplotypeCaller" + 
                    " -nct " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " --standard_min_confidence_threshold_for_calling 20" +
                    " --dbsnp " + WrapperUtility.ConvertWindowsPath(dbsnp) +
                    " -o " + WrapperUtility.ConvertWindowsPath(newVcf),
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

        private static string allGRCh37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz";
        private static string commonGRCh37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/common_all_20170710.vcf.gz";
        private static string allGRCh38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz";
        private static string commonGRCh38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz";

        public static void DownloadKnownSites(string binDirectory, string targetDirectory, bool commonOnly, bool GRCh37, bool GRCh38, out string knownSitesFilename)
        {
            knownSitesFilename = "";

            if (!GRCh37 && !GRCh38) return;

            string targetFile = commonOnly ?
                GRCh37 ? commonGRCh37 : commonGRCh38 :
                GRCh37 ? allGRCh37 : allGRCh38;

            knownSitesFilename = targetFile.Split('/').Last();
            knownSitesFilename = Path.GetFileNameWithoutExtension(knownSitesFilename);
            string newKnownSites = Path.GetFileNameWithoutExtension(knownSitesFilename) + ".ensembl.vcf";

            if (File.Exists(Path.Combine(targetDirectory, newKnownSites)))
            {
                knownSitesFilename = newKnownSites;
                return;
            }

            string scriptPath = Path.Combine(binDirectory, "scripts", "download_known_variants.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(targetDirectory),
                "wget " + targetFile,
                "gunzip " + knownSitesFilename + ".gz",
                "rm " + knownSitesFilename + ".gz"
            }).WaitForExit();
            File.Delete(scriptPath);

            Dictionary<string, string> chromMappings = File.ReadAllLines(GRCh37 ? Path.Combine(binDirectory, "ChromosomeMappings", "GRCh37_UCSC2ensembl.txt") : Path.Combine(binDirectory, "ChromosomeMappings", "GRCh38_UCSC2ensembl.txt"))
                .ToDictionary(line => line.Split('\t')[0], line => line.Split('\t')[1]);

            using (StreamReader reader = new StreamReader(Path.Combine(targetDirectory, knownSitesFilename)))
            using (StreamWriter writer = new StreamWriter(Path.Combine(targetDirectory, newKnownSites)))
            {
                while (true)
                {
                    string a = reader.ReadLine();
                    if (a == null) break;
                    string[] line = a.Split('\t');
                    if (chromMappings.TryGetValue(line[0], out string chr)) line[0] = chr;
                    writer.Write(String.Join("\t", line) + '\n');
                }
            }
            knownSitesFilename = newKnownSites;
            //File.Delete(Path.Combine(target_directory, known_sites_filename));
        }

        public static void Install(string currentDirectory)
        {
            if (Directory.Exists(Path.Combine(currentDirectory, "GenomeAnalysisTK.jar"))) return;
            string scriptPath = Path.Combine(currentDirectory, "scripts", "install_gatk.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                "while [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(currentDirectory, "GenomeAnalysisTK-*")) + " ]",
                "do",
                "  echo \"Genome Analysis Toolkit (GATK) not found.\nPlease download GATK from their website https://software.broadinstitute.org/gatk/download/ \nThen, place the file (.tar.bz2) in the folder " + 
                    currentDirectory + "\n\"",
                "  read -n 1 -s -r -p \"Press any key to continue\n\"",
                "done",
                "tar -jxvf GenomeAnalysisTK-*.tar.bz2",
                "rm GenomeAnalysisTK-*.tar.bz2",
                "mv GenomeAnalysisTK-*/GenomeAnalysisTK.jar .",
                "rm -r GenomeAnalysisTK-*",
                "git clone https://github.com/dpryan79/ChromosomeMappings.git",
            }).WaitForExit();
        }
    }
}
