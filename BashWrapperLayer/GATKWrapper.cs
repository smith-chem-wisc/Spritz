using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace ToolWrapperLayer
{
    public class GATKWrapper
    {

        #region Private Fields

        private static string allGRCh37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz";

        private static string commonGRCh37 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/common_all_20170710.vcf.gz";

        private static string allGRCh38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz";

        private static string commonGRCh38 = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz";

        private static Regex getFastaHeaderSequenceName = new Regex(@"(>)([\w\d\.\-]+)(.+)");

        //http://seqanswers.com/forums/showthread.php?t=22504
        private static string bamToChrBamOneLiner = " | awk 'BEGIN{FS=OFS=\"\t\"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3=\"chr\"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | ";

        #endregion Private Fields

        #region Public Fields

        public static string GATK = "java -Xmx20G -jar GenomeAnalysisTK.jar";
        public static string PICARD = "java -Xmx20G -jar picard.jar";

        #endregion Public Fields

        #region Public Methods

        public static void SubsetBam(string binDirectory, int threads, string bam, string genomeFasta, string genomeRegion, string outputBam)
        {
            string script_name = Path.Combine(binDirectory, "scripts", "subset_bam.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GATK +
                    " -T PrintReads" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(outputBam) +
                    " -L " + genomeRegion,
            }).WaitForExit();
        }

        /// <summary>
        /// Groups and sorts reads, and marks duplicates using Picard Tools. Then, splits and trims reads with SplitNCigarReads.
        /// 
        /// Run in parallel for each alignment, check that the RAM usage is okay... 
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="bam"></param>
        /// <param name="newBam"></param>
        public static void PrepareBamAndFasta(string binDirectory, int threads, string bam, string genomeFasta, string reference, out string newBam, out string ucscGenomeFasta)
        {
            // check if sorted and grouped and rename chromosomes
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            Dictionary<string, string> chromMappings = File.ReadAllLines(downloadGrch37 ?
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh37_ensembl2UCSC.txt") :
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh38_ensembl2UCSC.txt"))
                .Select(line => line.Split('\t'))
                .Where(x => x.Length > 1)
                .ToDictionary(line => line[0], line => line[1]);

            using (StreamWriter writer = new StreamWriter(Path.Combine(binDirectory, "scripts", "convertChromNamesOnTheFly.py")))
            {
                writer.Write("input sys\n");
                writer.Write("for line in sys.stdin\n");
                writer.Write("  new = line\n");
                foreach (var a in chromMappings.Where(x => x.Key.Length > 2 && x.Value.Length > 0)) // anything longer than MT is one of the weird contigs
                {
                    writer.Write("  new = new.replace(new, chr" + a.Key + ", " + a.Value + ")\n");
                }
                writer.Write("  sys.stdout.write(new)\n");
            }

            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerSorted");
            string readGroupedCheckfile = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerReadGrouped");
            newBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".UCSC.bam");
            string scriptPath = Path.Combine(binDirectory, "scripts", "check_sorted.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep '^@RG' > " + WrapperUtility.ConvertWindowsPath(readGroupedCheckfile),
                "samtools view -h " + WrapperUtility.ConvertWindowsPath(bam) + bamToChrBamOneLiner + "python scripts/convertChromNamesOnTheFly.py | samtools view -bS - > " + WrapperUtility.ConvertWindowsPath(newBam),
            }).WaitForExit();
            bool sorted = new FileInfo(Path.Combine(binDirectory, sortedCheckPath)).Length > 0;
            bool grouped = new FileInfo(Path.Combine(binDirectory, readGroupedCheckfile)).Length > 0;

            // reheader genomeFasta with UCSC chromsome names
            ucscGenomeFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".UCSC.fa");
            using (StreamReader reader = new StreamReader(Path.Combine(binDirectory, genomeFasta)))
            using (StreamWriter writer = new StreamWriter(Path.Combine(binDirectory, ucscGenomeFasta)))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null) break;
                    if (!line.StartsWith(">"))
                        writer.Write(line + '\n');
                    string sequenceName = getFastaHeaderSequenceName.Match(line).Groups[2].Value;
                    if (chromMappings.TryGetValue(sequenceName, out string chr))
                        writer.Write(getFastaHeaderSequenceName.Replace(line, m => m.Groups[1] + chr + m.Groups[3] + '\n'));
                }
            }

            // run commands for grouping or sorting
            string groupAndMaybeSortCommand;
            string groupSortBam = newBam;
            if (grouped && sorted) return;
            else if (!grouped)
            {
                groupSortBam = Path.Combine(Path.GetDirectoryName(newBam), Path.GetFileNameWithoutExtension(newBam) + (!sorted ? ".sorted.grouped.bam" : ".grouped.bam"));
                groupAndMaybeSortCommand = "if [ ! -f " + WrapperUtility.ConvertWindowsPath(groupSortBam) + " ]; then picard-tools AddOrReplaceReadGroups PU=platform  PL=illumina SM=sample LB=library" +
                    " I=" + WrapperUtility.ConvertWindowsPath(newBam) +
                    " O=" + WrapperUtility.ConvertWindowsPath(groupSortBam) +
                    (!sorted ? " SO=coordinate" : "") + 
                    "; fi";
            }
            else // grouped, but not sorted
            {
                groupSortBam = Path.Combine(Path.GetDirectoryName(newBam), Path.GetFileNameWithoutExtension(newBam) + ".sorted.bam");
                groupAndMaybeSortCommand = "if [ ! -f " + WrapperUtility.ConvertWindowsPath(groupSortBam) + " ]; then picard-tools SortSam SO=coordinate" +
                    " I=" + WrapperUtility.ConvertWindowsPath(newBam) +
                    " O=" + WrapperUtility.ConvertWindowsPath(groupSortBam) +
                    "; fi";
            }

            string markedDuplicatesBam = Path.Combine(Path.GetDirectoryName(groupSortBam), Path.GetFileNameWithoutExtension(groupSortBam) + ".marked.bam");
            string markedDuplicateMetrics = Path.Combine(Path.GetDirectoryName(groupSortBam), Path.GetFileNameWithoutExtension(groupSortBam) + ".marked.metrics");
            string splitTrimBam = Path.Combine(Path.GetDirectoryName(markedDuplicatesBam), Path.GetFileNameWithoutExtension(markedDuplicatesBam) + ".split.bam");
            string mapQReassigned = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".mapqfixed.bam");
            string scriptName2 = Path.Combine(binDirectory, "scripts", "picard.bash");

            string splitNCigarReadsCmd =
                    GATK +
                    " -T SplitNCigarReads" +
                    //" --num_threads " + threads.ToString() + // not supported
                    " -R " + WrapperUtility.ConvertWindowsPath(ucscGenomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " -U ALLOW_N_CIGAR_READS";

            WrapperUtility.GenerateAndRunScript(scriptName2, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                //reheaderCommand,
                groupAndMaybeSortCommand,
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ]]; then " +
                    "picard-tools MarkDuplicates" +
                    " I=" + WrapperUtility.ConvertWindowsPath(groupSortBam) +
                    " O=" + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) +
                    " M=" + WrapperUtility.ConvertWindowsPath(markedDuplicateMetrics) +
                    " AS=true; fi", // assume sorted
                "if [ -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ]; then rm " + WrapperUtility.ConvertWindowsPath(groupSortBam) + "; fi", // conserve space by removing former BAM

                GenomeFastaIndexCommand(ucscGenomeFasta),
                GenomeDictionaryIndexCommand(ucscGenomeFasta),

                "samtools index " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam),

                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ]]; then " + splitNCigarReadsCmd + " -fixMisencodedQuals; fi", // some datasets are probably going to have misencoded quality scores; this just subtracts 31 from all quality scores if possible...
                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ]]; then " + splitNCigarReadsCmd + "; fi",  // if it didn't run, it probably just found a correctly encoded read, so ditch the fixMisencodedQuals option

                "if [ -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ]; then " + 
                    "rm -f " +
                    WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " " +
                    WrapperUtility.ConvertWindowsPath(Path.Combine(Path.GetDirectoryName(markedDuplicatesBam), Path.GetFileNameWithoutExtension(markedDuplicatesBam) + ".bai")) +
                    "; fi",

                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(mapQReassigned) + " ] || [ ! -s " + WrapperUtility.ConvertWindowsPath(mapQReassigned) + " ]; then " +
                    GATK +
                    " -T PrintReads" +
                    " -R " + WrapperUtility.ConvertWindowsPath(ucscGenomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(mapQReassigned) +
                    " -rf ReassignMappingQuality; fi", // default mapping quality is 60; required for RNA-Seq aligners
                "samtools index " + WrapperUtility.ConvertWindowsPath(mapQReassigned),
                "if [ -f " + WrapperUtility.ConvertWindowsPath(mapQReassigned) + " ]; then rm " +
                    WrapperUtility.ConvertWindowsPath(splitTrimBam) + "; fi",
            }).WaitForExit();
            newBam = mapQReassigned;

            // run commands for marking duplicates and trimming reads
            File.Delete(sortedCheckPath);
            File.Delete(readGroupedCheckfile);
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
            string realignerTable = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".forIndelRealigner.intervals");
            newBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".realigned.bam");
            string script_name = Path.Combine(binDirectory, "scripts", "realign_indels.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(realignerTable) + " || ! -s " + WrapperUtility.ConvertWindowsPath(realignerTable) + " ]]; then " +
                    GATK +
                    " -T RealignerTargetCreator" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -o " +  WrapperUtility.ConvertWindowsPath(realignerTable) + 
                    "; fi",

                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(newBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(newBam) + " ]]; then " +
                    GATK +
                    " -T IndelRealigner" +
                    //" --num_threads " + threads.ToString() + // this tool can't do threaded analysis
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -targetIntervals " +  WrapperUtility.ConvertWindowsPath(realignerTable) +
                    " -o " + WrapperUtility.ConvertWindowsPath(newBam) + 
                    "; fi",
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
            string scriptPath = Path.Combine(binDirectory, "scripts", "base_recalibration.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(recalibrationTablePath) + " ]; then " +
                    GATK +
                    " -T BaseRecalibrator" +
                    //" --num_threads " + threads.ToString() + // doesn't support threaded runs
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -knownSites " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -o " + WrapperUtility.ConvertWindowsPath(recalibrationTablePath) + 
                    "; fi",
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
        public static void VariantCalling(string binDirectory, int threads, string genomeFasta, string bam, string dbsnpReferenceVcfPath, out string newVcf)
        {
            newVcf = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".vcf");
            string scriptName = Path.Combine(binDirectory, "scripts", "variant_calling.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(newVcf) + " || " + " ! -s " + WrapperUtility.ConvertWindowsPath(newVcf) + " ]]; then " +
                    GATK +
                    " -T HaplotypeCaller" +
                    " -nct " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " --standard_min_confidence_threshold_for_calling 20" +
                    " --dbsnp " + WrapperUtility.ConvertWindowsPath(dbsnpReferenceVcfPath) +
                    " -o " + WrapperUtility.ConvertWindowsPath(newVcf) + 
                    "; fi",
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

        public static void DownloadUCSCKnownVariantSites(string binDirectory, string targetDirectory, bool commonOnly, string reference, out string ucscKnownSitesPath)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            string targetFileLocation = commonOnly ?
                downloadGrch37 ? commonGRCh37 : commonGRCh38 :
                downloadGrch37 ? allGRCh37 : allGRCh38;
            string ucscKnownSitesFilename = targetFileLocation.Split('/').Last();
            ucscKnownSitesPath = Path.Combine(targetDirectory, Path.GetFileNameWithoutExtension(ucscKnownSitesFilename));


            if ((downloadGrch37 || downloadGrch38) 
                && !File.Exists(ucscKnownSitesPath))
            {
                string scriptPath = Path.Combine(binDirectory, "scripts", "downloadUcscVariants.bash");
                WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
                {
                    "cd " + WrapperUtility.ConvertWindowsPath(targetDirectory),
                    "wget " + targetFileLocation,
                    "gunzip " + WrapperUtility.ConvertWindowsPath(ucscKnownSitesPath) + ".gz",
                    "rm " +  WrapperUtility.ConvertWindowsPath(ucscKnownSitesPath) + ".gz"
                }).WaitForExit();
            }
        }

        public static void SortVCF(string binDirectory, string vcfPath, string genomeFastaPath, string sortedVcfPath)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "sortVcf.bash");
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath) + ".dict");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeDictionaryIndexCommand(genomeFastaPath),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(sortedVcfPath) + " ]; then " +
                    PICARD + " SortVcf I=" + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " O=" + WrapperUtility.ConvertWindowsPath(sortedVcfPath) +
                    " SEQUENCE_DICTIONARY=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) + 
                    "; fi"
            }).WaitForExit();
        }

        public static void Install(string currentDirectory)
        {
            string scriptPath = Path.Combine(currentDirectory, "scripts", "install_gatk.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(currentDirectory),
                "while [ ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(currentDirectory, "GenomeAnalysisTK*")) + " ]",
                "do",
                "  echo \"Genome Analysis Toolkit (GATK) not found.\nPlease download GATK from their website https://software.broadinstitute.org/gatk/download/ \nThen, place the file (.tar.bz2) in the folder " +
                    currentDirectory + "\n\"",
                "  read -n 1 -s -r -p \"Press any key to continue\n\"",
                "done",
                "if [ ! -f GenomeAnalysisTK.jar ]; then tar -jxvf GenomeAnalysisTK-*.tar.bz2; fi",
                "if [ ! -f GenomeAnalysisTK.jar ]; then rm GenomeAnalysisTK-*.tar.bz2; fi",
                "if [ ! -f GenomeAnalysisTK.jar ]; then mv GenomeAnalysisTK-*/GenomeAnalysisTK.jar .; fi",
                "if [ -f GenomeAnalysisTK.jar ]; then rm -r GenomeAnalysisTK-*; fi",
                "if [ ! -d ChromosomeMappings ]; then git clone https://github.com/dpryan79/ChromosomeMappings.git; fi",
                "if [ ! -f picard.jar ]; then wget https://github.com/broadinstitute/picard/releases/download/2.15.0/picard.jar; fi",
            }).WaitForExit();
        }

        #endregion Public Methods

        #region Private Methods

        private static string GenomeFastaIndexCommand(string genomeFastaPath)
        {
            return "samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFastaPath);
        }

        private static string GenomeDictionaryIndexCommand(string genomeFastaPath)
        {
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath) + ".dict");
            return "if [ -f " + WrapperUtility.ConvertWindowsPath(dictionaryPath) + " ]; then rm " + WrapperUtility.ConvertWindowsPath(dictionaryPath) + "; fi\n" +
                "picard-tools CreateSequenceDictionary R=" + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + " O=" + WrapperUtility.ConvertWindowsPath(dictionaryPath);
        }

        #endregion Private Methods

    }
}
