using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;

namespace ToolWrapperLayer
{
    public class STARWrapper
    {

        #region Public Properties

        public static string BamFileSuffix { get; } = "Aligned.out.bam";

        public static string SpliceJunctionFileSuffix { get; } = "SJ.out.tab";

        public static string LogFileSuffix { get; } = "Log.out";

        public static string LogFinalFileSuffix { get; } = "Log.final.out";

        public static string ChimericSamFileSuffix { get; } = "Chimeric.out.sam";

        public static string ChimericJunctionsFileSuffix { get; } = "Chimeric.out.junction";

        #endregion Public Properties

        public static void GenerateGenomeIndex(string binDirectory, int threads, string genomeDir, IEnumerable<string> genomeFastas, string geneModelGtfOrGff, int junctionOverhang = 100)
        {
            if (!Directory.Exists(genomeDir)) Directory.CreateDirectory(genomeDir);
            string fastas = String.Join(" ", genomeFastas.Select(f => WrapperUtility.ConvertWindowsPath(f)));
            string arguments =
                " --runMode genomeGenerate" +
                " --runThreadN " + threads.ToString() +
                " --genomeDir " + WrapperUtility.ConvertWindowsPath(genomeDir) + 
                " --genomeFastaFiles " + fastas +
                " --sjdbGTFfile " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGff) +
                (Path.GetExtension(geneModelGtfOrGff).StartsWith(".gff") ? " --sjdbGTFtagExonParentTranscript Parent" : "") +
                " --sjdbOverhang " + junctionOverhang.ToString();

            string script_name = Path.Combine(binDirectory, "scripts", "generate_genome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
        }

        public static void LoadGenome(string binDirectory, string genomeDir)
        {
            string arguments = " --genomeLoad " + STARGenomeLoadOption.LoadAndExit.ToString() +
                " --genomeDir '" + WrapperUtility.ConvertWindowsPath(genomeDir) + "'";
            string script_name = Path.Combine(binDirectory, "scripts", "load_genome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
        }


        public static void RemoveGenome(string binDirectory, string genomeDir)
        {
            string script_name = Path.Combine(binDirectory, "scripts", "removeGenome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() + " --genomeDir " + WrapperUtility.ConvertWindowsPath(genomeDir)
            }).WaitForExit();
        }

        // fastqs must have \n line endings, not \r\n
        public static void BasicAlignReads(string binDirectory, int threads, string genomeDir, string[] fastqFiles, string outprefix, bool strandSpecific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
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
                " --outSAMtype BAM Unsorted" +
                (strandSpecific ? "" : " --outSAMstrandField intronMotif") +
                " --chimSegmentMin 12" +
                " --chimJunctionOverhangMin 12" +
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) +
                read_command;

            string script_name = Path.Combine(binDirectory, "scripts", "align_reads.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + BamFileSuffix + " ]; then STAR/source/STAR" + arguments + "; fi",
                File.Exists(outprefix + BamFileSuffix) && genomeLoad == STARGenomeLoadOption.LoadAndRemove ? "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() : ""
            }).WaitForExit();
        }

        public void SuccessiveAlignReads(string binDirectory, int threads, string genomeDir, List<string[]> fastqFiles, string outprefix, bool strandSpecific = true)
        {
            LoadGenome(binDirectory, genomeDir);
            foreach(var fqs in fastqFiles)
            {
                BasicAlignReads(binDirectory, threads, genomeDir, fqs, outprefix, strandSpecific, STARGenomeLoadOption.LoadAndKeep);
            }
            RemoveGenome(binDirectory, genomeDir);
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

        public static void Install(string binDirectory)
        {
            bool downloadStar = !Directory.Exists(Path.Combine(binDirectory, "STAR"));
            string script_path = Path.Combine(binDirectory, "scripts", "installStar.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "git clone https://github.com/lh3/seqtk.git",
                "cd seqtk; make",
                "cd ..",
                downloadStar ? "git clone https://github.com/alexdobin/STAR.git" : "",
                downloadStar ? "cd STAR/source" : "",
                downloadStar ? "make STAR" : "",
                downloadStar ? "cd .." : "",
                downloadStar ? "export PATH=$PATH:" + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "STAR", "source")) : "",
            }).WaitForExit();
        }
    }
}
