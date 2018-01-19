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

        private static string allGRCh37UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz";

        private static string commonGRCh37UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/common_all_20170710.vcf.gz";

        // Malformed with empty alleles
        private static string GRCh37Ensembl = "ftp://ftp.ensembl.org/pub/release-75//variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz";

        private static string allGRCh38UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz";

        private static string commonGRCh38UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz";

        private static string GRCh38Ensembl = "ftp://ftp.ensembl.org/pub/release-81//variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz";

        private static Regex getFastaHeaderSequenceName = new Regex(@"(>)([\w\d\.\-]+)(.+)");

        private static Regex getISequenceHeaderSequenceName = new Regex(@"([\w\d\.\-]+)(.+)");

        //http://seqanswers.com/forums/showthread.php?t=22504
        private static string bamToChrBamOneLiner = " | awk 'BEGIN{FS=OFS=\"\t\"} (/^@/ && !/@SQ/){print $0} $2~/^SN:[1-9]|^SN:X|^SN:Y|^SN:MT/{print $0}  $3~/^[1-9]|X|Y|MT/{$3=\"chr\"$3; print $0} ' | sed 's/SN:/SN:chr/g' | sed 's/chrMT/chrM/g' | ";

        #endregion Private Fields

        #region General Public Methods

        public static string Gatk()
        {
            var performance = new PerformanceCounter("Memory", "Available MBytes");
            var memory = performance.NextValue();
            return "gatk/gatk --java-options -Xmx" + Math.Floor(memory) + "M";
        }

        #endregion General Public Methods

        #region Random Public Methods

        public static void SubsetBam(string binDirectory, int threads, string bam, string genomeFasta, string genomeRegion, string outputBam)
        {
            string script_name = Path.Combine(binDirectory, "scripts", "subset_bam.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                Gatk() +
                    " PrintReads" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(outputBam) +
                    " -L " + genomeRegion,
            }).WaitForExit();
        }

        public static void SortVCF(string binDirectory, string vcfPath, string genomeFastaPath, string sortedVcfPath)
        {
            string tmpDir = Path.Combine(binDirectory, "tmp");
            Directory.CreateDirectory(tmpDir);
            string scriptPath = Path.Combine(binDirectory, "scripts", "sortVcf.bash");
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath) + ".dict");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeDictionaryIndexCommand(genomeFastaPath),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(sortedVcfPath) + " ]; then " +
                    Gatk() + // formerly picard
                    " SortVcf I=" + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " O=" + WrapperUtility.ConvertWindowsPath(sortedVcfPath) +
                    " SEQUENCE_DICTIONARY=" + WrapperUtility.ConvertWindowsPath(dictionaryPath) +
                    " TMP_DIR=" + WrapperUtility.ConvertWindowsPath(tmpDir) +
                    "; fi"
            }).WaitForExit();
        }

        #endregion Random Public Methods

        #region Downloading Reference Variant Databases

        public static void DownloadUCSCKnownVariantSites(string binDirectory, string targetDirectory, bool commonOnly, string reference, out string ucscKnownSitesPath)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            string targetFileLocation =
                commonOnly ?
                    (downloadGrch37 ? commonGRCh37UCSC : commonGRCh38UCSC) :
                    (downloadGrch37 ? allGRCh37UCSC : allGRCh38UCSC);
            string ucscKnownSitesFilename = targetFileLocation.Split('/').Last();
            ucscKnownSitesPath = Path.Combine(targetDirectory, Path.GetFileNameWithoutExtension(ucscKnownSitesFilename));

            if ((downloadGrch37 || downloadGrch38) && !File.Exists(ucscKnownSitesPath))
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

        public static void DownloadEnsemblKnownVariantSites(string binDirectory, string targetDirectory, bool commonOnly, string reference, out string ensemblKnownSitesPath)
        {
            DownloadUCSCKnownVariantSites(binDirectory, targetDirectory, commonOnly, reference, out string ucscKnownSitesPath);
            ConvertVCFChromosomesUCSC2Ensembl(binDirectory, ucscKnownSitesPath, reference, out ensemblKnownSitesPath);
        }

        public static void ConvertVCFChromosomesUCSC2Ensembl(string binDirectory, string vcfPath, string reference, out string ensemblVcfPath)
        {
            ensemblVcfPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath) + ".ensembl.vcf");
            if (File.Exists(ensemblVcfPath)) return;
            Dictionary<string, string> ucsc2EnsemblChromosomeMappings = EnsemblDownloadsWrapper.UCSC2EnsemblChromosomeMappings(binDirectory, reference);
            using (StreamReader reader = new StreamReader(vcfPath))
            using (StreamWriter writer = new StreamWriter(ensemblVcfPath))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null)
                    {
                        break;
                    }
                    if (line.StartsWith("#"))
                    {
                        writer.Write(line + "\n");
                    }
                    string[] splitLine = line.Split('\t');
                    if (splitLine.Length > 1 && ucsc2EnsemblChromosomeMappings.TryGetValue(splitLine[0], out string newChrom))
                    {
                        splitLine[0] = newChrom;
                        writer.Write(String.Join("\t", splitLine) + '\n');
                    }
                }
            }
        }

        #endregion Downloading Reference Variant Databases

        #region Variant Calling

        /// <summary>
        /// HaplotypeCaller for calling variants on each RNA-Seq BAM file individually
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="dedupedBam"></param>
        /// <param name=""></param>
        public static void VariantCalling(string binDirectory, int threads, string genomeFasta, string dedupedBam, string dbsnpReferenceVcfPath, out string newVcf)
        {
            string fixedQualsBam = Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam) + ".fixedQuals.bam");
            string splitTrimBam = Path.Combine(Path.GetDirectoryName(fixedQualsBam), Path.GetFileNameWithoutExtension(fixedQualsBam) + ".split.bam");
            string unfliteredVcf = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".vcf");
            newVcf = Path.Combine(Path.GetDirectoryName(unfliteredVcf), Path.GetFileNameWithoutExtension(unfliteredVcf) + "filtered.vcf");

            string fixMisencodedQualsCmd =
                Gatk() +
                " FixMisencodedBaseQualityReads" +
                " -I " + WrapperUtility.ConvertWindowsPath(dedupedBam) +
                " -O " + WrapperUtility.ConvertWindowsPath(fixedQualsBam);

            string splitNCigarReadsCmd1 =
                Gatk() +
                " SplitNCigarReads" +
                //" --num_threads " + threads.ToString() + // not supported
                " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                " -I " + WrapperUtility.ConvertWindowsPath(fixedQualsBam) +
                " -O " + WrapperUtility.ConvertWindowsPath(splitTrimBam)
                //" -rf ReassignOneMappingQuality" + // doing this with STAR
                //" -RMQF 255" +
                //" -RMQT 60" + // default mapping quality is 60; required for RNA-Seq aligners
                //" -U ALLOW_N_CIGAR_READS"
                ;

            string splitNCigarReadsCmd2 =
                Gatk() +
                " SplitNCigarReads" +
                //" --num_threads " + threads.ToString() + // not supported
                " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                " -I " + WrapperUtility.ConvertWindowsPath(dedupedBam) +
                " -O " + WrapperUtility.ConvertWindowsPath(splitTrimBam)
                //" -rf ReassignOneMappingQuality" + // doing this with STAR
                //" -RMQF 255" +
                //" -RMQT 60" + // default mapping quality is 60; required for RNA-Seq aligners
                //" -U ALLOW_N_CIGAR_READS"
                ;

            string scriptName = Path.Combine(binDirectory, "scripts", "variantCalling.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                // split and trim reads (some datasets are probably going to have misencoded quality scores; -fixMisencodedQuals just subtracts 31 from all quality scores if possible...)
                // if it doesn't run with -fixMisencodedQuals, it probably just found a correctly encoded read, so ditch the fixMisencodedQuals option
                "samtools index " + WrapperUtility.ConvertWindowsPath(dedupedBam),
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(fixedQualsBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(fixedQualsBam) + " ) ]]; then " + fixMisencodedQualsCmd + "; fi",
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ) ]]; then " + splitNCigarReadsCmd1 + "; fi",
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ) ]]; then " + splitNCigarReadsCmd2 + "; fi",
                "samtools index " + WrapperUtility.ConvertWindowsPath(splitTrimBam),

                // call variants
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(unfliteredVcf) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(unfliteredVcf) + " ]; then " +
                    Gatk() +
                    " HaplotypeCaller" +
                    " --native-pair-hmm-threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " --min-base-quality-score 20" +
                    " --dont-use-soft-clipped-bases true" + // for RNA-Seq
                    " --dbsnp " + WrapperUtility.ConvertWindowsPath(dbsnpReferenceVcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(unfliteredVcf) +
                    "; fi",

                // filter variants (RNA-Seq specific params... need to check out recommendations before using DNA-Seq)
                //"if [ ! -f " + WrapperUtility.ConvertWindowsPath(newVcf) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(newVcf) + " ]; then " +
                //    Gatk() +
                //    " -T VariantFiltration" +
                //    " -nct " + threads.ToString() +
                //    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                //    " -V " + WrapperUtility.ConvertWindowsPath(unfliteredVcf) +
                //    " -window 35 -cluster 3" + // filter out clusters of 3 snps within 35 bases (https://software.broadinstitute.org/gatk/documentation/topic?name=methods)
                //    " -filterName FS -filter \"FS > 30.0\"" +
                //    " -filterName QD -filter \"QD < 2.0\"" +
                //    " -o " + WrapperUtility.ConvertWindowsPath(newVcf) +
                //    "; fi",
            }).WaitForExit();
            newVcf = unfliteredVcf;
        }

        #endregion Variant Calling

        #region Defunct Methods, or might be used for DNA-Seq

        /// <summary>
        /// Groups and sorts reads, and marks duplicates using Picard Tools. Then, splits and trims reads with SplitNCigarReads.
        /// 
        /// TODO: separate out the ucsc conversions from before
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="bam"></param>
        /// <param name="newBam"></param>
        public static void PrepareBamAndFasta(string binDirectory, int threads, string bam, string genomeFasta, string reference, out string newBam, bool convertToUCSC = true)
        {
            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerSorted");
            string readGroupedCheckfile = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerReadGrouped");
            string sortedBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".sorted.bam");
            string groupedBam = Path.Combine(Path.GetDirectoryName(sortedBam), Path.GetFileNameWithoutExtension(sortedBam) + ".grouped.bam");
            string markedDuplicatesBam = Path.Combine(Path.GetDirectoryName(groupedBam), Path.GetFileNameWithoutExtension(groupedBam) + ".marked.bam");
            string markedDuplicateMetrics = Path.Combine(Path.GetDirectoryName(groupedBam), Path.GetFileNameWithoutExtension(groupedBam) + ".marked.metrics");

            string tmpDir = Path.Combine(binDirectory, "tmp");
            Directory.CreateDirectory(tmpDir);
            string scriptName2 = Path.Combine(binDirectory, "scripts", "picard." + Path.GetFileNameWithoutExtension(bam) + ".bash");
            WrapperUtility.GenerateAndRunScript(scriptName2, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),

                GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep SO:coordinate > " + WrapperUtility.ConvertWindowsPath(sortedCheckPath),
                "samtools view -H " + WrapperUtility.ConvertWindowsPath(bam) + " | grep '^@RG' > " + WrapperUtility.ConvertWindowsPath(readGroupedCheckfile),
                
                // group and sort (note, using picard-tools works, but picard.jar somehow is trucating the BAM files)
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(groupedBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(groupedBam) + " ) && " +
                    " ( ! -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ) ]]; then " +
                    Gatk() +
                    " AddOrReplaceReadGroups" + 
                    " -PU platform  -PL illumina -SM sample -LB library" +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " -O " + WrapperUtility.ConvertWindowsPath(groupedBam) +
                    " -SO coordinate" +
                    " --TMP_DIR " + WrapperUtility.ConvertWindowsPath(tmpDir) +
                    "; fi",
                "if [[ -f " + WrapperUtility.ConvertWindowsPath(groupedBam) + " && -s " + WrapperUtility.ConvertWindowsPath(groupedBam) + " ]]; then rm " + WrapperUtility.ConvertWindowsPath(sortedBam) + "; fi", // conserve space by removing former BAM

                // mark duplicates (AS means assume sorted; note, using picard-tools works, but picard.jar somehow is trucating the BAM files)
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ) ]]; then " +
                    Gatk() + " MarkDuplicates" + // formerly picard
                    " -I " + WrapperUtility.ConvertWindowsPath(groupedBam) +
                    " -O " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) +
                    " -M " + WrapperUtility.ConvertWindowsPath(markedDuplicateMetrics) +
                    " --TMP_DIR " + WrapperUtility.ConvertWindowsPath(tmpDir) +
                    " -AS true" +
                    "; fi",
                "if [[ -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " && -s " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ]]; then rm " + WrapperUtility.ConvertWindowsPath(groupedBam) + "; fi", // conserve space by removing former BAM
                "samtools index " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam),

            }).WaitForExit();
            newBam = markedDuplicatesBam;

            // run commands for marking duplicates and trimming reads
            File.Delete(sortedCheckPath);
            File.Delete(readGroupedCheckfile);
        }

        /// <summary>
        /// Realigns indels for a given BAM file
        /// 
        /// This is no longer a required step for HaploytypeCaller, used for variant calling
        /// 
        /// Needs updating after switching to GATK 4.0 ... no longer called IndelRealigner, and it looks like the targetIntervals is gone
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
                    Gatk() +
                    " RealignerTargetCreator" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -o " +  WrapperUtility.ConvertWindowsPath(realignerTable) +
                    "; fi",

                "if [[ ! -f " + WrapperUtility.ConvertWindowsPath(newBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(newBam) + " ]]; then " +
                    Gatk() +
                    " LeftAlignIndels" +
                    //" --num_threads " + threads.ToString() + // this tool can't do threaded analysis
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -known " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -targetIntervals " +  WrapperUtility.ConvertWindowsPath(realignerTable) +
                    " --OUTPUT " + WrapperUtility.ConvertWindowsPath(newBam) +
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
                    Gatk() +
                    " BaseRecalibrator" +
                    //" --num_threads " + threads.ToString() + // doesn't support threaded runs
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -knownSites " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -O " + WrapperUtility.ConvertWindowsPath(recalibrationTablePath) +
                    "; fi",
            }).WaitForExit();
        }

        // using GenotypeGVCFs on all VCF files together, for DNA sequence analysis
        public static Process genotype()
        {
            return null;
        }

        // using VariantRecalibrator, ApplyRecalibration, and then VariantFiltration -- apparently mostly for DNAseq
        public static Process filter_variants()
        {
            return null;
        }

        #endregion Defunct Methods, or might be used for DNA-Seq

        #region Private Methods

        private static string GenomeFastaIndexCommand(string genomeFastaPath)
        {
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + ".fai ]; then samtools faidx " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + "; fi";
        }

        private static string GenomeDictionaryIndexCommand(string genomeFastaPath)
        {
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath) + ".dict");
            return "if [ ! -f " + WrapperUtility.ConvertWindowsPath(dictionaryPath) + " ]; then " + //rm " + WrapperUtility.ConvertWindowsPath(dictionaryPath) + "; fi\n" +
                Gatk() + // formerly picard
                    " CreateSequenceDictionary" + 
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + 
                    " -O " + WrapperUtility.ConvertWindowsPath(dictionaryPath) +
                    "; fi";
        }

        #endregion Private Methods

    }
}

// from PrepareBam originally
// check if sorted and grouped and rename chromosomes
//newBam = bam;
//string ucscBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".UCSC.bam");
//if (convertToUCSC && !File.Exists(ucscBam))
//{
//    newBam = ucscBam;
//    Dictionary<string, string> chromMappings = EnsemblDownloadsWrapper.Ensembl2UCSCChromosomeMappings(binDirectory, reference);

//    // future speedup: could check BAM header to only convert the chroms that are there
//    // future speedup: divide up amongst ProcessorCount # of python scripts and pipe them together
//    using (StreamWriter writer = new StreamWriter(Path.Combine(binDirectory, "scripts", "convertChromNamesOnTheFly.py")))
//    {
//        writer.Write("import sys\n");
//        writer.Write("for line in sys.stdin:\n");
//        writer.Write("  new = line\n");
//        foreach (var a in chromMappings.Where(x => x.Key.Length > 2 && x.Value.Length > 0)) // anything longer than MT is one of the weird contigs
//        {
//            writer.Write("  new = new.replace(\"chr" + a.Key + "\", \"" + a.Value + "\")\n");
//        }
//        writer.Write("  sys.stdout.write(new)\n");
//    }

//    // reheader genomeFasta with UCSC chromsome names; remove the ones that aren't in the UCSC universe
//    Genome ucscGenome = new Genome(genomeFasta);
//    ucscGenome.Chromosomes = ucscGenome.Chromosomes
//        .Where(x => chromMappings.TryGetValue(getISequenceHeaderSequenceName.Match(x.ID).Groups[1].Value, out string chr) && chr.Length > 0).ToList();
//    foreach (ISequence chrom in ucscGenome.Chromosomes)
//    {
//        string sequenceName = getISequenceHeaderSequenceName.Match(chrom.ID).Groups[1].Value;
//        if (chromMappings.TryGetValue(sequenceName, out string chr))
//        {
//            chrom.ID = getISequenceHeaderSequenceName.Replace(chrom.ID, m => chr + m.Groups[2]);
//        }
//    }

//    ucscGenomeFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".UCSC.fa");
//    Genome.WriteFasta(ucscGenome.KaryotypicOrder(), ucscGenomeFasta);
//}

//// sort
//"if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(sortedBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(sortedBam) + " ) && " +
//    " ( ! -f " + WrapperUtility.ConvertWindowsPath(groupedBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(groupedBam) + " ) && " +
//    " ( ! -f " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(markedDuplicatesBam) + " ) && " +
//    " ( ! -f " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(splitTrimBam) + " ) && " +
//    " ( ! -f " + WrapperUtility.ConvertWindowsPath(mapQReassigned) + " || ! -s " + WrapperUtility.ConvertWindowsPath(mapQReassigned) + " ) ]]; then " +
//    "samtools sort -f -@ " + Environment.ProcessorCount.ToString() + " -m " + Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue()) + "M " +
//        WrapperUtility.ConvertWindowsPath(bam) + " " + WrapperUtility.ConvertWindowsPath(sortedBam) +
//        "; fi",