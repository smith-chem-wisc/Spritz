using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;

namespace RNASeqAnalysisWrappers
{
    public class STARWrapper
    {
        public static void generate_genome_index(string bin_directory, int threads, string genomeDir, IEnumerable<string> genomeFastas, string gene_model, int junctionOverhang = 100)
        {
            if (!Directory.Exists(genomeDir)) Directory.CreateDirectory(genomeDir);
            string fastas = String.Join(" ", genomeFastas.Select(f => WrapperUtility.convert_windows_path(f)));
            string arguments =
                " --runMode genomeGenerate" +
                " --runThreadN " + threads.ToString() +
                " --genomeDir " + WrapperUtility.convert_windows_path(genomeDir) + 
                " --genomeFastaFiles " + fastas +
                " --sjdbGTFfile " + WrapperUtility.convert_windows_path(gene_model) +
                (Path.GetExtension(gene_model).StartsWith(".gff") ? " --sjdbGTFtagExonParentTranscript Parent" : "") +
                " --sjdbOverhang " + junctionOverhang.ToString();

            string script_name = Path.Combine(bin_directory, "generate_genome.sh");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
            File.Delete(script_name);
        }

        public static void load_genome(string bin_directory, string genomeDir)
        {
            string arguments = " --genomeLoad LoadAndExit" +
                " --genomeDir '" + WrapperUtility.convert_windows_path(genomeDir) + "'";
            string script_name = Path.Combine(bin_directory, "load_genome.sh");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
            File.Delete(script_name);
        }


        public static void remove_genome(string bin_directory)
        {
            string script_name = Path.Combine(bin_directory, "remove_genome.sh");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR --genomeLoad Remove"
            }).WaitForExit();
            File.Delete(script_name);
        }

        public static void basic_align_reads(string bin_directory, int threads, string genomeDir, string[] fastq_files, string outprefix, string genomeLoad = "NoSharedMemory", bool strand_specific = true)
        {
            string reads_in = "\"" + String.Join("\" \"", fastq_files.Select(f => WrapperUtility.convert_windows_path(f))) + "\"";
            string read_command = fastq_files.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastq_files.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --genomeLoad " + genomeLoad +
                " --genomeDir \"" + WrapperUtility.convert_windows_path(genomeDir) + "\"" +
                " --readFilesIn " + reads_in +
                " --outSAMtype BAM Unsorted" +
                " --outFileNamePrefix " + WrapperUtility.convert_windows_path(outprefix) +
                read_command;

            string script_name = Path.Combine(bin_directory, "align_reads.sh");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
            //File.Delete(script_name);
        }

        public void successive_align_reads(string bin_directory, int threads, string genomeDir, List<string[]> fastq_files, string outprefix, bool strand_specific = true)
        {
            load_genome(bin_directory, genomeDir);
            foreach(var fqs in fastq_files)
            {
                basic_align_reads(bin_directory, threads, genomeDir, fqs, outprefix, "LoadAndKeep", strand_specific);
            }
            remove_genome(bin_directory);
        }

        public static bool check_strand_specificity(string bin_directory, int threads, string genomeDir, string[] fastq_files, string bed_genes, string current_directory)
        {
            List<string> new_files = new List<string>();
            foreach(string file in fastq_files)
            {
                string new_path = Path.Combine(current_directory, Path.GetFileNameWithoutExtension(file) + ".segment.fastq");
                new_files.Add(new_path);
                using (StreamWriter writer = new StreamWriter(new_path))
                {
                    using (FileStream fstream = new FileStream(file, FileMode.Open, FileAccess.Read, FileShare.Read))
                    {
                        Stream stream = Path.GetExtension(file) == ".gz" ?
                            (Stream)(new GZipStream(fstream, CompressionMode.Decompress)) :
                             fstream;
                        StreamReader reader = new StreamReader(stream);

                        for (int i = 0; i < 4000000 && reader.Peek() > -1; i++)
                        {
                            writer.WriteLine(reader.ReadLine());
                        }
                    }
                }
            }

            basic_align_reads(bin_directory, threads, genomeDir, new_files.ToArray(), 
                Path.Combine(Path.GetDirectoryName(new_files[0]), Path.GetFileNameWithoutExtension(new_files[0]) + "."));

            bool result = RSeQCWrapper.check_strand_specificity(bin_directory, Path.Combine(Path.GetDirectoryName(new_files[0]), Path.GetFileNameWithoutExtension(new_files[0]) + ".Aligned.out.bam"), bed_genes, current_directory);
            foreach (string file in Directory.GetFiles(current_directory, "*.segment.fastq").Concat(Directory.GetFiles(current_directory, "*strand_specificity.*")))
            {
                File.Delete(file);
            }
            return result;
        }

        public static void install(string current_directory)
        {
            if (Directory.Exists(Path.Combine(current_directory, "STAR"))) return;
            string script_path = Path.Combine(current_directory, "install_star.sh");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "git clone https://github.com/alexdobin/STAR.git",
                "cd STAR/source",
                "make STAR",
                "git submodule update --init --recursive", // include STAR-Fusion
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
