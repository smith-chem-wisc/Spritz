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


        public static void GenerateGenomeIndex(string bin_directory, int threads, string genomeDir, IEnumerable<string> genomeFastas, string geneModelGtfOrGff, int junctionOverhang = 100)
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

            string script_name = Path.Combine(bin_directory, "scripts", "generate_genome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
        }

        public static void LoadGenome(string bin_directory, string genomeDir)
        {
            string arguments = " --genomeLoad " + STARGenomeLoadOption.LoadAndExit.ToString() +
                " --genomeDir '" + WrapperUtility.ConvertWindowsPath(genomeDir) + "'";
            string script_name = Path.Combine(bin_directory, "scripts", "load_genome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
        }


        public static void RemoveGenome(string bin_directory)
        {
            string script_name = Path.Combine(bin_directory, "scripts", "remove_genome.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString()
            }).WaitForExit();
        }

        // fastqs must have \n line endings, not \r\n
        public static void BasicAlignReads(string bin_directory, int threads, string genomeDir, string[] fastq_files, string outprefix, bool strand_specific = true, STARGenomeLoadOption genomeLoad = STARGenomeLoadOption.NoSharedMemory)
        {
            string reads_in = "\"" + String.Join("\" \"", fastq_files.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\"";
            string read_command = fastq_files.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastq_files.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --genomeLoad " + genomeLoad.ToString() +
                " --runThreadN " + threads.ToString() +
                " --genomeDir \"" + WrapperUtility.ConvertWindowsPath(genomeDir) + "\"" +
                " --readFilesIn " + reads_in +
                " --outSAMtype BAM Unsorted" +
                (strand_specific ? "" : " --outSAMstrandField intronMotif") +
                " --chimSegmentMin 12" +
                " --chimJunctionOverhangMin 12" +
                " --outFileNamePrefix " + WrapperUtility.ConvertWindowsPath(outprefix) +
                read_command;

            string script_name = Path.Combine(bin_directory, "scripts", "align_reads.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "if [ ! -f " + WrapperUtility.ConvertWindowsPath(outprefix) + BamFileSuffix + " ]; then STAR/source/STAR" + arguments + "; fi",
                File.Exists(outprefix + BamFileSuffix) && genomeLoad == STARGenomeLoadOption.LoadAndRemove ? "STAR/source/STAR --genomeLoad " + STARGenomeLoadOption.Remove.ToString() : ""
            }).WaitForExit();
        }

        public void SuccessiveAlignReads(string bin_directory, int threads, string genomeDir, List<string[]> fastq_files, string outprefix, bool strand_specific = true)
        {
            LoadGenome(bin_directory, genomeDir);
            foreach(var fqs in fastq_files)
            {
                BasicAlignReads(bin_directory, threads, genomeDir, fqs, outprefix, strand_specific, STARGenomeLoadOption.LoadAndKeep);
            }
            RemoveGenome(bin_directory);
        }

        public static void SubsetFastqs(string[] fastq_files, int reads, string current_directory, out string[] new_files)
        {
            // Note: fastqs must have \n line endings, not \r\n

            List<string> new_ones = new List<string>();
            foreach(string file in fastq_files)
            {
                string new_path = Path.Combine(Path.GetDirectoryName(file), Path.GetFileNameWithoutExtension(file) + ".segment.fastq");
                new_ones.Add(new_path);

                using (StreamWriter writer = new StreamWriter(new_path))
                using (FileStream fstream = new FileStream(file, FileMode.Open, FileAccess.Read, FileShare.Read))
                {
                    Stream stream = Path.GetExtension(file) == ".gz" ?
                        (Stream)(new GZipStream(fstream, CompressionMode.Decompress)) :
                            fstream;
                    StreamReader reader = new StreamReader(stream);

                    for (int i = 0; i < reads * 4 && reader.Peek() > -1; i++)
                    {
                        writer.Write(reader.ReadLine() + '\n');
                    }
                }
            }
            new_files = new_ones.ToArray();
        }

        public static void Install(string binDirectory)
        {
            bool downloadStar = !Directory.Exists(Path.Combine(binDirectory, "STAR"));
            string script_path = Path.Combine(binDirectory, "scripts", "installStar.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                downloadStar ? "git clone https://github.com/alexdobin/STAR.git" : "",
                downloadStar ? "cd STAR/source" : "",
                downloadStar ? "make STAR" : "",
                downloadStar ? "cd .." : "",
                downloadStar ? "export PATH=$PATH:" + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "STAR", "source")) : "",
            }).WaitForExit();
        }
    }
}
