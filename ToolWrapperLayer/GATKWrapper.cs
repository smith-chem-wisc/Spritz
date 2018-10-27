using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace ToolWrapperLayer
{
    /// <summary>
    /// GATK is the genome analysis toolkit, used for calling variants from nucleotide sequencing data and for preparing alignment maps for variant calling.
    /// </summary>
    public class GATKWrapper :
        IInstallable
    {
        #region dbSNP URLs

        /// <summary>
        /// All dbSNP reference alleles for GRCh37, used for variant calling
        /// </summary>
        private static string AllGRCh37UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/All_20170710.vcf.gz";

        /// <summary>
        /// Common dbSNP reference alleles for GRCh37, used for variant calling
        /// </summary>
        private static string CommonGRCh37UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/GATK/common_all_20170710.vcf.gz";

        /// <summary>
        /// dbSNP reference alleles for GRCh37 from Ensembl. Don't use this. The VCF is malformed with some empty alleles.
        /// </summary>
        private static string GRCh37Ensembl = "ftp://ftp.ensembl.org/pub/release-75//variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz";

        /// <summary>
        /// All dbSNP reference alleles for GRCh38, used for variant calling
        /// </summary>
        private static string AllGRCh38UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/All_20170710.vcf.gz";

        /// <summary>
        /// Common dbSNP reference alleles for GRCh38, used for variant calling
        /// </summary>
        private static string CommonGRCh38UCSC = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/GATK/common_all_20170710.vcf.gz";

        /// <summary>
        /// dbSNP reference alleles for GRCh37 from Ensembl. Don't use this. The VCF might be malformed with some empty alleles, although I haven't checked like I did for the GRCh37 one.
        /// </summary>
        private static string GRCh38Ensembl = "ftp://ftp.ensembl.org/pub/release-81//variation/vcf/homo_sapiens/Homo_sapiens.vcf.gz";

        #endregion dbSNP URLs

        private static string Version = "4.0.11.0";
        private static Regex getFastaHeaderSequenceName = new Regex(@"(>)([\w\d\.\-]+)(.+)");
        private static Regex getISequenceHeaderSequenceName = new Regex(@"([\w\d\.\-]+)(.+)");
        private string TargetFileLocation;
        public string UcscKnownSitesPath { get; private set; }
        public string EnsemblKnownSitesPath { get; private set; }
        public string SplitTrimBamPath { get; private set; }
        public string HaplotypeCallerGvcfPath { get; private set; } // intermediate file with all genomic loci
        public string HaplotypeCallerVcfPath { get; private set; }
        public string FilteredHaplotypeCallerVcfPath { get; private set; }
        public string PreparedBamPath { get; private set; }
        public string RecalibrationTablePath { get; private set; }
        public string RecalibratedBamPath { get; private set; }
        public string SortedVcfPath { get; private set; }

        /// <summary>
        /// Writes an installation script for installing GATK from pre-built binaries.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallGatk.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d gatk ]; then",
                $"  wget --no-check https://github.com/broadinstitute/gatk/releases/download/{Version}/gatk-{Version}.zip",
                $"  unzip gatk-{Version}.zip",
                $"  rm gatk-{Version}.zip",
                $"  mv gatk-{Version} gatk",
                "fi",
                "if [ ! -d ChromosomeMappings ]; then git clone https://github.com/dpryan79/ChromosomeMappings.git; fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing gatk.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        /// <summary>
        /// Generic method for subsetting a BAM file. Useful for testing new methods.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="bam"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="genomeRegion"></param>
        /// <param name="outputBam"></param>
        public void SubsetBam(string spritzDirectory, string analysisDirectory, int threads, string bam, string genomeFasta, string genomeRegion, string outputBam)
        {
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "SubsetBam.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                Gatk() +
                    " PrintReads" +
                    " --num_threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " -o " + WrapperUtility.ConvertWindowsPath(outputBam) +
                    " -L " + genomeRegion,
            }).WaitForExit();
        }

        /// <summary>
        /// Sorts a VCF file based on the chromosome ordering in the genome fasta.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="vcfPath"></param>
        /// <param name="genomeFastaPath"></param>
        /// <param name="sortedVcfPath"></param>
        public List<string> SortVCF(string spritzDirectory, string analysisDirectory, string vcfPath, string genomeFastaPath)
        {
            SortedVcfPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath)) + ".sorted.vcf";
            string tmpDir = Path.Combine(spritzDirectory, "tmp");
            Directory.CreateDirectory(tmpDir);
            string dictionaryPath = Path.Combine(Path.GetDirectoryName(genomeFastaPath), Path.GetFileNameWithoutExtension(genomeFastaPath) + ".dict");
            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                GenomeDictionaryIndexCommand(genomeFastaPath),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(SortedVcfPath) + " ]; then " +
                    Gatk() + // formerly picard
                    " SortVcf -I " + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(SortedVcfPath) +
                    " --SEQUENCE_DICTIONARY " + WrapperUtility.ConvertWindowsPath(dictionaryPath) +
                    " --TMP_DIR " + WrapperUtility.ConvertWindowsPath(tmpDir) +
                    "; fi"
            };
        }

        /// <summary>
        /// Reference dbSNP vcf file exists
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="commonOnly"></param>
        /// <param name="reference"></param>
        /// <returns></returns>
        public bool KnownVariantSitesFileExists(string spritzDirectory, bool commonOnly, string reference)
        {
            bool downloadGrch37 = string.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool downloadGrch38 = string.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            TargetFileLocation =
                commonOnly ?
                    (downloadGrch37 ? CommonGRCh37UCSC : CommonGRCh38UCSC) :
                    (downloadGrch37 ? AllGRCh37UCSC : AllGRCh38UCSC);
            string ucscKnownSitesFilename = TargetFileLocation.Split('/').Last();
            UcscKnownSitesPath = Path.Combine(spritzDirectory, Path.GetFileNameWithoutExtension(ucscKnownSitesFilename));
            return (downloadGrch37 || downloadGrch38) && File.Exists(UcscKnownSitesPath);
        }

        /// <summary>
        /// Downloads dbSNP reference VCF file if it doesn't exist
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="commonOnly"></param>
        /// <param name="reference"></param>
        /// <returns></returns>
        public string DownloadUCSCKnownVariantSites(string spritzDirectory, bool commonOnly, string reference)
        {
            if (!KnownVariantSitesFileExists(spritzDirectory, commonOnly, reference))
            {
                WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(spritzDirectory, "DownloadUcscVariants.bash"), new List<string>
                {
                    "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                    "wget " + TargetFileLocation,
                    "gunzip " + WrapperUtility.ConvertWindowsPath(UcscKnownSitesPath) + ".gz",
                    "rm " +  WrapperUtility.ConvertWindowsPath(UcscKnownSitesPath) + ".gz"
                }).WaitForExit();
            }
            return UcscKnownSitesPath;
        }

        public string DownloadEnsemblKnownVariantSites(string spritzDirectory, bool commonOnly, string reference)
        {
            DownloadUCSCKnownVariantSites(spritzDirectory, commonOnly, reference);
            EnsemblKnownSitesPath = ConvertVCFChromosomesUCSC2Ensembl(spritzDirectory, UcscKnownSitesPath, reference);

            // indexing is used for most GATK tools
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(spritzDirectory, "IndexKnownVariantSites.bash"), new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(UcscKnownSitesPath) + ".idx ]; then " + Gatk() + " IndexFeatureFile -F " + WrapperUtility.ConvertWindowsPath(UcscKnownSitesPath) + "; fi",
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(EnsemblKnownSitesPath) + ".idx ]; then " + Gatk() + " IndexFeatureFile -F " + WrapperUtility.ConvertWindowsPath(EnsemblKnownSitesPath) + "; fi",
            }).WaitForExit();
            return EnsemblKnownSitesPath;
        }

        public string ConvertVCFChromosomesUCSC2Ensembl(string spritzDirectory, string vcfPath, string reference)
        {
            string convertedEnsemblVcfPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath) + ".ensembl.vcf");
            if (File.Exists(convertedEnsemblVcfPath)) { return convertedEnsemblVcfPath; }
            Dictionary<string, string> ucsc2EnsemblChromosomeMappings = EnsemblDownloadsWrapper.UCSC2EnsemblChromosomeMappings(spritzDirectory, reference);
            using (StreamReader reader = new StreamReader(vcfPath))
            using (StreamWriter writer = new StreamWriter(convertedEnsemblVcfPath))
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
                        writer.Write(string.Join("\t", splitLine) + '\n');
                    }
                }
            }
            return convertedEnsemblVcfPath;
        }

        /// <summary>
        /// Splits and trims reads splice junction reads with SplitNCigarReads.
        /// Apparently cigars are genomic intervals, and splice junctions are represented by a bunch of N's (unkonwn nucleotide), HaplotypeCaller requires splitting them in the BAM file.
        ///
        /// It's tempting to want to run a few of these at the same time because it's not well parallelized. It's just not worth it. It uses quite a bit of RAM and racks the I/O at the beginning when reading the BAM files.
        /// Could possibly do 4 at a time on 128 GB RAM and 28 processors.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="dedupedBam"></param>
        /// <param name="splitTrimBam"></param>
        /// <returns></returns>
        public List<string> SplitNCigarReads(string spritzDirectory, string genomeFasta, string dedupedBam)
        {
            string fixedQualsBam = Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam) + ".fixedQuals.bam");
            SplitTrimBamPath = Path.Combine(Path.GetDirectoryName(fixedQualsBam), Path.GetFileNameWithoutExtension(fixedQualsBam) + ".split.bam");

            // This also filters malformed reads
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
                " -O " + WrapperUtility.ConvertWindowsPath(SplitTrimBamPath)
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
                " -O " + WrapperUtility.ConvertWindowsPath(SplitTrimBamPath)
                //" -rf ReassignOneMappingQuality" + // doing this with STAR
                //" -RMQF 255" +
                //" -RMQT 60" + // default mapping quality is 60; required for RNA-Seq aligners
                //" -U ALLOW_N_CIGAR_READS"
                ;

            List<string> commands = new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                SamtoolsWrapper.GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                // split and trim reads (some datasets are probably going to have misencoded quality scores; -fixMisencodedQuals just subtracts 31 from all quality scores if possible...)
                // exit code of 2 means that the FixMisencodedQualityBaseReads errored out because there were correctly encode base quality scores
                SamtoolsWrapper.IndexBamCommand(dedupedBam),
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(fixedQualsBam) + " || ! -s " + WrapperUtility.ConvertWindowsPath(fixedQualsBam) + " ) ]]; then",
                "  " + fixMisencodedQualsCmd,
                "  if [ $? -ne 2 ]; then",
                "    " + splitNCigarReadsCmd1,
                "  else",
                "    " + splitNCigarReadsCmd2,
                "  fi",
                "fi",
                SamtoolsWrapper.IndexBamCommand(SplitTrimBamPath),
            };
            return commands;
        }

        /// <summary>
        /// HaplotypeCaller for calling variants on each RNA-Seq BAM file individually.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="splitTrimBam"></param>
        /// <param name="dbsnpReferenceVcfPath"></param>
        /// <param name="newVcf"></param>
        public List<string> VariantCalling(string spritzDirectory, ExperimentType experimentType, int threads, string genomeFasta, string splitTrimBam, string dbsnpReferenceVcfPath)
        {
            HaplotypeCallerGvcfPath = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".g.vcf.gz");
            HaplotypeCallerVcfPath = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".g.gt.vcf");
            FilteredHaplotypeCallerVcfPath = Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + ".g.gt.NoIndels.vcf");
            var vcftools = new VcfToolsWrapper();

            List<string> commands = new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                SamtoolsWrapper.GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                // check that reference VCF is indexed
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(dbsnpReferenceVcfPath) + ".idx ]; then " + Gatk() + " IndexFeatureFile -F " + WrapperUtility.ConvertWindowsPath(dbsnpReferenceVcfPath) + "; fi",

                // call variants
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) + " ]; then " +
                    Gatk() +
                    " HaplotypeCaller" +
                    " --native-pair-hmm-threads " + threads.ToString() +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(splitTrimBam) +
                    " --min-base-quality-score 20" +
                    (experimentType == ExperimentType.RNASequencing ? " --dont-use-soft-clipped-bases true" : "") +
                    " --dbsnp " + WrapperUtility.ConvertWindowsPath(dbsnpReferenceVcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) +
                    " -ERC GVCF" + // this prompts phasing!
                    " --max-mnp-distance 3" + // note: this can't be used for joint genotyping here, but this setting is available in mutect2 for doing tumor vs normal calls
                "; fi",

                // index compressed gvcf file
                $"if [ ! -f {WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath)}.idx ]; then {Gatk()} IndexFeatureFile -F {WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath)}; fi",

                // genotype the gvcf file into a traditional vcf file
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ]; then " +
                    Gatk() +
                    " GenotypeGVCFs" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -V " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) +
                "; fi",

                // filter out indels
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ]; then " +
                    Gatk() +
                    " SelectVariants" +
                    " --select-type-to-exclude INDEL" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -V " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(FilteredHaplotypeCallerVcfPath) +
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
            };
            return commands;
        }

        public List<string> CombineAndGenotypeGvcfs(string spritzDirectory, string genomeFasta, List<string> gvcfPaths)
        {
            if (gvcfPaths == null || gvcfPaths.Count <= 1)
            {
                throw new ArgumentException("CombineAndGenotypeGvcfs exception: no gvcfs were specified to combine");
            }
            int uniqueSuffix = 1;
            foreach (string f in gvcfPaths)
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            HaplotypeCallerGvcfPath = Path.Combine(Path.GetDirectoryName(gvcfPaths.First()), $"combined{uniqueSuffix}.g.vcf.gz");
            HaplotypeCallerVcfPath = Path.Combine(Path.GetDirectoryName(HaplotypeCallerGvcfPath), $"{Path.GetFileNameWithoutExtension(Path.GetFileNameWithoutExtension(HaplotypeCallerGvcfPath))}.gt.vcf");
            FilteredHaplotypeCallerVcfPath = Path.Combine(Path.GetDirectoryName(HaplotypeCallerGvcfPath), $"{Path.GetFileNameWithoutExtension(HaplotypeCallerGvcfPath)}.NoIndels.vcf");

            List<string> commands = new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                SamtoolsWrapper.GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta)
            };

            foreach (string gvcf in gvcfPaths)
            {
                // double check that the compressed gvcf file is indexed
                commands.Add($"if [ ! -f {WrapperUtility.ConvertWindowsPath(gvcf)}.idx ]; then {Gatk()} IndexFeatureFile -F {WrapperUtility.ConvertWindowsPath(gvcf)}; fi");
            }

            // combine GVCFs
            string combineCommand =
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) + " ]; then " +
                    Gatk() +
                    " CombineGVCFs" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -V " + string.Join(" -V ", gvcfPaths.Select(gvcf => WrapperUtility.ConvertWindowsPath(gvcf))) +
                    " -O " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) +
                "; fi";
            commands.Add(combineCommand);
            commands.Add($"if [ ! -f {WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath)}.idx ]; then {Gatk()} IndexFeatureFile -F {WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath)}; fi");

            // genotype the gvcf file into a traditional vcf file
            string genotypeCommand =
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ]; then " +
                    Gatk() +
                    " GenotypeGVCFs" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -V " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerGvcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) +
                "; fi";
            commands.Add(genotypeCommand);

            // filter out indels
            string filterIndelsCommand =
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ] || [ " + " ! -s " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) + " ]; then " +
                    Gatk() +
                    " SelectVariants" +
                    " --select-type-to-exclude INDEL" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -V " + WrapperUtility.ConvertWindowsPath(HaplotypeCallerVcfPath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(FilteredHaplotypeCallerVcfPath) +
                "; fi";
            commands.Add(filterIndelsCommand);

            return commands;
        }

        #region Defunct Methods, or might be used for DNA-Seq

        /// <summary>
        /// Groups (I'm just using one group, so it's more a formality) and sorts reads. Marks duplicates.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="bam"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="reference"></param>
        /// <param name="newBam"></param>
        /// <param name="convertToUCSC"></param>
        public List<string> PrepareBamAndFasta(string spritzDirectory, string analysisDirectory, int threads, string bam, string genomeFasta, string reference)
        {
            string sortedCheckPath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerSorted");
            string readGroupedCheckfile = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".headerReadGrouped");
            string sortedBam = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".sorted.bam");
            string groupedBam = Path.Combine(Path.GetDirectoryName(sortedBam), Path.GetFileNameWithoutExtension(sortedBam) + ".grouped.bam");
            string markedDuplicatesBam = Path.Combine(Path.GetDirectoryName(groupedBam), Path.GetFileNameWithoutExtension(groupedBam) + ".marked.bam");
            string markedDuplicateMetrics = Path.Combine(Path.GetDirectoryName(groupedBam), Path.GetFileNameWithoutExtension(groupedBam) + ".marked.metrics");

            string tmpDir = Path.Combine(spritzDirectory, "tmp");
            Directory.CreateDirectory(tmpDir);
            string scriptName2 = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "picard." + Path.GetFileNameWithoutExtension(bam) + ".bash");
            List<string> commands = new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),

                SamtoolsWrapper.GenomeFastaIndexCommand(genomeFasta),
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
            };
            PreparedBamPath = markedDuplicatesBam;

            // run commands for marking duplicates and trimming reads
            File.Delete(sortedCheckPath);
            File.Delete(readGroupedCheckfile);
            return commands;
        }

        /// <summary>
        /// Creates recalibration table and recalibrates reads.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="genomeFasta"></param>
        /// <param name="bam"></param>
        /// <param name="recalibrationTablePath"></param>
        /// <param name="knownSitesVcf"></param>
        public List<string> BaseRecalibration(string spritzDirectory, string analysisDirectory, string genomeFasta, string bam, string knownSitesVcf)
        {
            RecalibrationTablePath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".recaltable");
            RecalibratedBamPath = Path.Combine(Path.GetDirectoryName(bam), Path.GetFileNameWithoutExtension(bam) + ".recal.bam");

            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                SamtoolsWrapper.GenomeFastaIndexCommand(genomeFasta),
                GenomeDictionaryIndexCommand(genomeFasta),

                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(RecalibrationTablePath) + " ]; then " +
                    Gatk() +
                    " BaseRecalibrator" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    (knownSitesVcf != "" ? " -knownSites " + WrapperUtility.ConvertWindowsPath(knownSitesVcf) : "") +
                    " -O " + WrapperUtility.ConvertWindowsPath(RecalibrationTablePath) +
                    "; fi",

                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(RecalibrationTablePath) + " ]; then " +
                    Gatk() +
                    " ApplyBQSR" +
                    " -R " + WrapperUtility.ConvertWindowsPath(genomeFasta) +
                    " -I " + WrapperUtility.ConvertWindowsPath(bam) +
                    " --bqsr-recal-file " + WrapperUtility.ConvertWindowsPath(RecalibrationTablePath) +
                    " -O " + WrapperUtility.ConvertWindowsPath(RecalibratedBamPath) +
                    "; fi",
            };
        }

        /// <summary>
        /// Generic command for calling GATK, allowing use of all free memory.
        /// </summary>
        /// <returns></returns>
        public static string Gatk()
        {
            var performance = new PerformanceCounter("Memory", "Available MBytes");
            var memory = performance.NextValue();
            return "gatk/gatk --java-options -Xmx" + Math.Floor(memory) + "M";
        }

        /// <summary>
        /// Writes an installation script for installing GATK from source. Requires root permissions to install gradle the first time.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static string WriteGitCloneInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallGatk.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -f gatk/build/libs/gatk.jar ]; then",
                "  git clone https://github.com/broadinstitute/gatk.git",
                "  cd gatk",
                "  ./gradlew localJar",
                "  cd ..",
                "fi",
                "if [ ! -d ChromosomeMappings ]; then git clone https://github.com/dpryan79/ChromosomeMappings.git; fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// using GenotypeGVCFs on all VCF files together, for DNA sequence analysis
        /// </summary>
        /// <returns></returns>
        public static Process genotype()
        {
            return null;
        }

        /// <summary>
        /// using VariantRecalibrator, ApplyRecalibration, and then VariantFiltration -- apparently mostly for DNAseq
        /// </summary>
        /// <returns></returns>
        public static Process filter_variants()
        {
            return null;
        }

        #endregion Defunct Methods, or might be used for DNA-Seq

        #region Private Methods

        /// <summary>
        /// Creates a dictionary for the genome fasta, used by many GATK tools.
        /// </summary>
        /// <param name="genomeFastaPath"></param>
        /// <returns></returns>
        private string GenomeDictionaryIndexCommand(string genomeFastaPath)
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