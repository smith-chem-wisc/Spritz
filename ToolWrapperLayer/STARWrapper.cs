using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// STAR is a fast and accurate spliced aligner for RNA-Seq data. It requires a lot of RAM, ~40 GB of free RAM.
    /// It has an option for two-pass alignment, which improves the accuracy of splice junction detection.
    /// </summary>
    public class STARWrapper :
        IInstallable
    {
        #region Public Properties

        /// <summary>
        /// Output BAM file. Suffix tagged onto the output prefix
        /// </summary>
        public static string BamFileSuffix { get; } = "Aligned.out.bam";

        /// <summary>
        /// Output sorted BAM file. Suffix tagged onto the output prefix
        /// </summary>
        public static string SortedBamFileSuffix { get; } = "Aligned.sortedByCoord.out.bam";

        /// <summary>
        /// Output deduped BAM file. Suffix tagged onto the output prefix
        /// </summary>
        public static string DedupedBamFileSuffix { get; } = "Aligned.sortedByCoord.outProcessed.out.bam";

        /// <summary>
        /// Output deduping log file. Suffix tagged onto the output prefix
        /// </summary>
        public static string DedupedBamFileLog { get; } = "Aligned.sortedByCoord.outLog.out";

        /// <summary>
        /// Output splice junction file. Suffix tagged onto the output prefix
        /// </summary>
        public static string SpliceJunctionFileSuffix { get; } = "SJ.out.tab";

        /// <summary>
        /// Log file. Suffix tagged onto the output prefix
        /// </summary>
        public static string LogFileSuffix { get; } = "Log.out";

        /// <summary>
        /// Final log file. Suffix tagged onto the output prefix
        /// </summary>
        public static string LogFinalFileSuffix { get; } = "Log.final.out";

        /// <summary>
        /// Chimeric alignment file. Suffix tagged onto the output prefix
        /// </summary>
        public static string ChimericSamFileSuffix { get; } = "Chimeric.out.sam";

        /// <summary>
        /// Chimeric junction file. Suffix tagged onto the output prefix
        /// </summary>
        public static string ChimericJunctionsFileSuffix { get; } = "Chimeric.out.junction";

        /// <summary>
        /// File name used to check for too many open files bug, followed by an integer
        /// </summary>
        public static string ThreadCheckFileSuffix { get; } = "Chimeric.out.junction.thread";

        #endregion Public Properties

        #region Genome Index Methods

        /// <summary>
        /// Generates indices for STAR alignments.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeDir"></param>
        /// <param name="genomeFastas"></param>
        /// <param name="geneModelGtfOrGff"></param>
        /// <param name="sjdbFileChrStartEnd"></param>
        /// <param name="junctionOverhang"></param>
        /// <returns></returns>
        public static List<string> GenerateGenomeIndex(string binDirectory, int threads, string genomeDir, IEnumerable<string> genomeFastas, string geneModelGtfOrGff, string sjdbFileChrStartEnd = "", int junctionOverhang = 100)
        {
            string fastas = String.Join(" ", genomeFastas.Select(f => WrapperUtility.ConvertWindowsPath(f)));
            string arguments =
                " --runMode genomeGenerate" +
                " --runThreadN " + threads.ToString() +
                " --genomeDir " + WrapperUtility.ConvertWindowsPath(genomeDir) +
                " --genomeFastaFiles " + fastas +
                " --sjdbGTFfile " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGff) +
                (File.Exists(sjdbFileChrStartEnd) ? " --limitSjdbInsertNsj 1200000 --sjdbFileChrStartEnd " + WrapperUtility.ConvertWindowsPath(sjdbFileChrStartEnd) : "") +
                (Path.GetExtension(geneModelGtfOrGff).StartsWith(".gff") ? " --sjdbGTFtagExonParentTranscript Parent" : "") +
                " --sjdbOverhang " + junctionOverhang.ToString();

            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -d " + WrapperUtility.ConvertWindowsPath(genomeDir) + " ]; then mkdir " + WrapperUtility.ConvertWindowsPath(genomeDir) + "; fi",
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(genomeDir, "SA")) + " || ! -s " + WrapperUtility.ConvertWindowsPath(Path.Combine(genomeDir, "SA")) + " ) ]]; then STAR/source/STAR" + arguments + "; fi"
            };
        }

        /// <summary>
        /// Removes genome indices from memory.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="genomeDir"></param>
        /// <returns></returns>
        public static List<string> RemoveGenome(string binDirectory, string genomeDir)
        {
            string script_name = Path.Combine(binDirectory, "scripts", "removeGenome.bash");
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() + " --genomeDir " + WrapperUtility.ConvertWindowsPath(genomeDir)
            };
        }

        #endregion Genome Index Methods

        #region Installation Methods

        /// <summary>
        /// Writes an installation script for STAR. Also installs seqtk, which is useful for subsetting fastq files.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installScripts", "installStar.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "git clone https://github.com/lh3/seqtk.git",
                "cd seqtk; make",
                "cd ..",
                "if [ ! -d STAR ]; then ",
                "  git clone https://github.com/alexdobin/STAR.git",
                "  cd STAR/source",
                "  make STAR",
                "  cd ..",
                "  export PATH=$PATH:" + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "STAR", "source")),
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing STAR and seqtk.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

        #region First-Pass Alignment Methods

        /// <summary>
        /// Aligns reads and outputs junctions of spliced alignments. Does not output an alignment map.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeDir"></param>
        /// <param name="fastqFiles"></param>
        /// <param name="outprefix"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="genomeLoad"></param>
        /// <returns></returns>
        public static List<string> FirstPassAlignmentCommands(string binDirectory, int threads, string genomeDir, string[] fastqFiles, string outprefix, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
        {
            return BasicAlignReadCommands(binDirectory, threads, genomeDir, fastqFiles, outprefix, strandSpecific, genomeLoad, "None");
        }

        /// <summary>
        /// Bundles splice jucntions from first pass alignments into a single splice junction file for second-pass alignment.
        /// Excludes splice junctions for mitochondrial chromosome alignments.
        /// </summary>
        /// <param name="spliceJunctionOuts"></param>
        /// <param name="uniqueSuffix"></param>
        /// <param name="spliceJunctionStarts"></param>
        /// <returns></returns>
        public static List<string> ProcessFirstPassSpliceCommands(List<string> spliceJunctionOuts, int uniqueSuffix, out string spliceJunctionStarts)
        {
            if (spliceJunctionOuts.Count == 0) throw new ArgumentException("STARWrapper.ProcessFirstPassSpliceCommands: No splice junctions detected for second-pass genome generation.");
            spliceJunctionStarts = Path.Combine(Path.GetDirectoryName(spliceJunctionOuts[0]), "combined" + uniqueSuffix.ToString() + "." + SpliceJunctionFileSuffix);
            return new List<string>
            {
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(spliceJunctionStarts) + " ]; then " +
                    "awk 'BEGIN {OFS=\"\t\"; strChar[0]=\".\"; strChar[1]=\"+\"; strChar[2]=\"-\";} {if($5>0){print $1,$2,$3,strChar[$4]}}' " +
                    String.Join(" ", spliceJunctionOuts.Select(f => WrapperUtility.ConvertWindowsPath(f))) +
                    " | grep -v 'MT' >> " +
                    WrapperUtility.ConvertWindowsPath(spliceJunctionStarts) +
                    "; fi"
            };
        }

        #endregion First-Pass Alignment Methods

        #region Public Methods

        /// <summary>
        /// Aligns reads and outputs alignment map and chimeric alignments.
        /// Note: fastqs must have \n line endings, not \r\n.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeDir"></param>
        /// <param name="fastqFiles"></param>
        /// <param name="outprefix"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="genomeLoad"></param>
        /// <param name="outSamType"></param>
        /// <returns></returns>
        public static List<string> BasicAlignReadCommands(string binDirectory, int threads, string genomeDir, string[] fastqFiles, string outprefix, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory, string outSamType = "BAM Unsorted")
        {
            string reads_in = "\"" + String.Join("\" \"", fastqFiles.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\"";
            string read_command = fastqFiles.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastqFiles.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --genomeLoad " + genomeLoad.ToString() +
                " --runThreadN " + threads.ToString() +
                " --genomeDir \"" + WrapperUtility.ConvertWindowsPath(genomeDir) + "\"" +
                " --readFilesIn " + reads_in +
                " --outSAMtype " + outSamType +
                " --limitBAMsortRAM " + (Math.Round(Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue() * 1e6), 0)).ToString() +
                (strandSpecific ? "" : " --outSAMstrandField intronMotif") + // puts in XS attribute for unstranded data, which is used by Cufflinks
                " --chimSegmentMin 12" +
                " --chimJunctionOverhangMin 12" +
                " --outFilterIntronMotifs RemoveNoncanonical" + // for cufflinks
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) +
                read_command;

            string fileToCheck = WrapperUtility.ConvertWindowsPath(outprefix) + (outSamType.Contains("Sorted") ? SortedBamFileSuffix : outSamType.Contains("Unsorted") ? BamFileSuffix : SpliceJunctionFileSuffix);
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(fileToCheck) + " || ! -s " + WrapperUtility.ConvertWindowsPath(fileToCheck) + " ) ]]; then STAR/source/STAR" + arguments + "; fi",
                File.Exists(outprefix + BamFileSuffix) && genomeLoad == STARGenomeLoadOption.LoadAndRemove ? "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() : ""
            };
        }

        /// <summary>
        /// Aligns reads and outputs alignment map and chimeric alignments. Duplicate reads are removed (deduped) from the alignment map, a step that's recommended for variant calling.
        /// Note: fastqs must have \n line endings, not \r\n.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="genomeDir"></param>
        /// <param name="fastqFiles"></param>
        /// <param name="outprefix"></param>
        /// <param name="strandSpecific"></param>
        /// <param name="genomeLoad"></param>
        /// <returns></returns>
        public static List<string> AlignRNASeqReadsForVariantCalling(string binDirectory, int threads, string genomeDir, string[] fastqFiles,
            string outprefix, bool overwriteStarAlignment, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
        {
            string reads_in = "\"" + String.Join("\" \"", fastqFiles.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\"";
            string read_command = fastqFiles.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastqFiles.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";

            string alignmentArguments =
                " --genomeLoad " + genomeLoad.ToString() +
                " --runMode alignReads" +
                " --runThreadN " + threads.ToString() +
                " --genomeDir \"" + WrapperUtility.ConvertWindowsPath(genomeDir) + "\"" +
                " --readFilesIn " + reads_in +
                " --outSAMtype BAM SortedByCoordinate" +
                " --limitBAMsortRAM " + (Math.Round(Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue() * 1e6), 0)).ToString() +
                (strandSpecific ? "" : " --outSAMstrandField intronMotif") + // puts in XS attribute for unstranded data, which is used by Cufflinks
                " --chimSegmentMin 12" +
                " --chimJunctionOverhangMin 12" +
                " --outFilterIntronMotifs RemoveNoncanonical" + // for cufflinks
                " --outSAMattrRGline ID:1 PU:platform  PL:illumina SM:sample LB:library" + // this could shorten the time for samples that aren't multiplexed in preprocessing for GATK
                " --outSAMmapqUnique 60" + // this could be used to ensure compatibility with GATK without having to use the GATK hacks
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) +
                read_command; // note in the future, two sets of reads can be comma separated here, and the RGline can also be comma separated to distinguish them later

            string dedupArguments =
                " --runMode inputAlignmentsFromBAM" +
                " --bamRemoveDuplicatesType UniqueIdentical" + // this could shorten the time for samples that aren't multiplexed, too; might only work with sortedBAM input from --inputBAMfile
                " --limitBAMsortRAM " + (Math.Round(Math.Floor(new PerformanceCounter("Memory", "Available MBytes").NextValue() * 1e6), 0)).ToString() +
                " --runThreadN " + threads.ToString() +
                " --inputBAMfile " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix +
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) + Path.GetFileNameWithoutExtension(SortedBamFileSuffix);

            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),

                overwriteStarAlignment ? "" : "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix + " || ! -s " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix + " ) ]]; then",
                "  STAR/source/STAR" + alignmentArguments,
                overwriteStarAlignment ? "" : "fi",
                SamtoolsWrapper.IndexBamCommand(binDirectory, WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix),

                overwriteStarAlignment ? "" : "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + DedupedBamFileSuffix + " || ! -s " + WrapperUtility.ConvertWindowsPath(outprefix) + DedupedBamFileSuffix + " ) ]]; then",
                "  STAR/source/STAR" + dedupArguments,
                overwriteStarAlignment ? "" : "fi",
                SamtoolsWrapper.IndexBamCommand(binDirectory, WrapperUtility.ConvertWindowsPath(outprefix) + DedupedBamFileSuffix),

                File.Exists(outprefix + BamFileSuffix) && File.Exists(outprefix + DedupedBamFileSuffix) && genomeLoad == STARGenomeLoadOption.LoadAndRemove ?
                    "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() :
                    "",
            };
        }

        /// <summary>
        /// Uses seqtk to get a subset of reads from a (pair of) fastq file(s).
        /// Note: fastqs must have \n line endings, not \r\n.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="fastqFiles"></param>
        /// <param name="reads"></param>
        /// <param name="currentDirectory"></param>
        /// <param name="newFfiles"></param>
        /// <param name="useSeed"></param>
        /// <param name="seed"></param>
        public static void SubsetFastqs(string binDirectory, string[] fastqFiles, int reads, string currentDirectory, out string[] newFfiles, bool useSeed = false, int seed = 0)
        {
            newFfiles = new string[] { Path.Combine(Path.GetDirectoryName(fastqFiles[0]), Path.GetFileNameWithoutExtension(fastqFiles[0]) + ".segment.fastq") };
            if (fastqFiles.Length > 1)
                newFfiles = new string[] { newFfiles[0], Path.Combine(Path.GetDirectoryName(fastqFiles[1]), Path.GetFileNameWithoutExtension(fastqFiles[1]) + ".segment.fastq") };

            string script_path = Path.Combine(binDirectory, "scripts", "subsetReads.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -s " + WrapperUtility.ConvertWindowsPath(newFfiles[0]) + " ]; then",
                "  echo \"Subsetting " + reads.ToString() + " reads from " + String.Join(",", fastqFiles) + "\"",
                "  seqtk/seqtk sample" + (useSeed || fastqFiles.Length > 1 ? " -s" + seed.ToString() : "") + " " + WrapperUtility.ConvertWindowsPath(fastqFiles[0]) + " " + reads.ToString() + " > " + WrapperUtility.ConvertWindowsPath(newFfiles[0]),
                fastqFiles.Length > 1 ? "  seqtk/seqtk sample -s" + seed.ToString() + " " + WrapperUtility.ConvertWindowsPath(fastqFiles[1]) + " " + reads.ToString() + " > " + WrapperUtility.ConvertWindowsPath(newFfiles[1]) : "",
                "fi"
            }).WaitForExit();
        }

        /// <summary>
        /// Gets the Windows-formatted path to the STAR executable
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public static string GetStarDirectoryPath(string binDirectory)
        {
            return Path.Combine(binDirectory, "STAR", "source");
        }

        #endregion Public Methods
    }
}