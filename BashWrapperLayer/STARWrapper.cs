using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Diagnostics;

namespace ToolWrapperLayer
{
    public class STARWrapper
    {

        #region Public Properties

        public static string BamFileSuffix { get; } = "Aligned.out.bam";

        public static string SortedBamFileSuffix { get; } = "Aligned.sortedByCoord.out.bam";

        public static string DedupedBamFileSuffix { get; } = "Aligned.sortedByCoord.outProcessed.out.bam";

        public static string DedupedBamFileLog { get; } = "Aligned.sortedByCoord.outLog.out";

        public static string SpliceJunctionFileSuffix { get; } = "SJ.out.tab";

        public static string LogFileSuffix { get; } = "Log.out";

        public static string LogFinalFileSuffix { get; } = "Log.final.out";

        public static string ChimericSamFileSuffix { get; } = "Chimeric.out.sam";

        public static string ChimericJunctionsFileSuffix { get; } = "Chimeric.out.junction";

        #endregion Public Properties

        #region Genome Index Methods

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
                "mkdir " + WrapperUtility.ConvertWindowsPath(genomeDir),
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(Path.Combine(genomeDir, "SA")) + " || ! -s " + WrapperUtility.ConvertWindowsPath(Path.Combine(genomeDir, "SA")) + " ) ]]; then STAR/source/STAR" + arguments + "; fi"
            };
        }

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

        #region First-Pass Alignment Methods

        public static List<string> FirstPassAlignmentCommands(string binDirectory, int threads, string genomeDir, string[] fastqFiles, string outprefix, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
        {
            return BasicAlignReadCommands(binDirectory, threads, genomeDir, fastqFiles, outprefix, strandSpecific, genomeLoad, "None");
        }

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

        #endregion

        // fastqs must have \n line endings, not \r\n
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

        // fastqs must have \n line endings, not \r\n
        public static List<string> AlignRNASeqReadsForVariantCalling(string binDirectory, int threads, string genomeDir, string[] fastqFiles, string outprefix, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
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
                " --runThreadN " + threads.ToString() +
                " --inputBAMfile " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix +
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) + Path.GetFileNameWithoutExtension(SortedBamFileSuffix);

            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix + " || ! -s " + WrapperUtility.ConvertWindowsPath(outprefix) + SortedBamFileSuffix + " ) ]]; then STAR/source/STAR" + alignmentArguments + "; fi",
                "if [[ ( ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + DedupedBamFileSuffix + " || ! -s " + WrapperUtility.ConvertWindowsPath(outprefix) + DedupedBamFileSuffix + " ) ]]; then STAR/source/STAR" + dedupArguments + "; fi",
                File.Exists(outprefix + BamFileSuffix) && File.Exists(outprefix + DedupedBamFileSuffix) && genomeLoad == STARGenomeLoadOption.LoadAndRemove ? "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() : ""
            };
        }

        public static void SubsetFastqs(string binDirectory, string[] fastqFiles, int reads, string currentDirectory, out string[] newFfiles, bool useSeed = false, int seed = 0)
        {
            // Note: fastqs must have \n line endings, not \r\n
            newFfiles = new string[] { Path.Combine(Path.GetDirectoryName(fastqFiles[0]), Path.GetFileNameWithoutExtension(fastqFiles[0]) + ".segment.fastq") };
            if (fastqFiles.Length > 1)
                newFfiles = new string[] { newFfiles[0], Path.Combine(Path.GetDirectoryName(fastqFiles[1]), Path.GetFileNameWithoutExtension(fastqFiles[1]) + ".segment.fastq") };

            string script_path = Path.Combine(binDirectory, "scripts", "subsetReads.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "seqtk/seqtk sample" + (useSeed || fastqFiles.Length > 1 ? " -s" + seed.ToString() : "") + " " + WrapperUtility.ConvertWindowsPath(fastqFiles[0]) + " " + reads.ToString() + " > " + WrapperUtility.ConvertWindowsPath(newFfiles[0]),
                fastqFiles.Length > 1 ? "seqtk/seqtk sample -s" + seed.ToString() + " " + WrapperUtility.ConvertWindowsPath(fastqFiles[1]) + " " + reads.ToString() + " > " + WrapperUtility.ConvertWindowsPath(newFfiles[1]) : "",
            }).WaitForExit();
        }

        public static string WriteInstallScript(string binDirectory)
        {
            string scriptPath = Path.Combine(binDirectory, "scripts", "installStar.bash");
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
    }
}
