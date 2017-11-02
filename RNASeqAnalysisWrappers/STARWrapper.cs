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

            string script_name = Path.Combine(bin_directory, "generate_genome.bash");
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
            string script_name = Path.Combine(bin_directory, "load_genome.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
            File.Delete(script_name);
        }


        public static void remove_genome(string bin_directory)
        {
            string script_name = Path.Combine(bin_directory, "remove_genome.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR --genomeLoad Remove"
            }).WaitForExit();
            File.Delete(script_name);
        }

        // fastqs must have \n line endings, not \r\n
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
                " --runThreadN " + threads.ToString() +
                " --genomeDir \"" + WrapperUtility.convert_windows_path(genomeDir) + "\"" +
                " --readFilesIn " + reads_in +
                " --outSAMtype BAM Unsorted" +
                (strand_specific ? "" : " --outSAMstrandField intronMotif") +
                " --chimSegmentMin 12" +
                " --chimJunctionOverhangMin 12" +
                " --outFileNamePrefix " + WrapperUtility.convert_windows_path(outprefix) +
                read_command;

            string script_name = Path.Combine(bin_directory, "align_reads.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR" + arguments
            }).WaitForExit();
            File.Delete(script_name);
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

        public static void subset_fastqs(string[] fastq_files, int reads, string current_directory, out string[] new_files)
        {
            List<string> new_ones = new List<string>();
            foreach(string file in fastq_files)
            {
                string new_path = Path.Combine(current_directory, Path.GetFileNameWithoutExtension(file) + ".segment.fastq");
                new_ones.Add(new_path);
                using (StreamWriter writer = new StreamWriter(new_path))
                {
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
            }
            new_files = new_ones.ToArray();
        }

        // fastqs must have \n line endings, not \r\n
        public static void star_fusion(string bin_directory, bool GRCh37, bool GRCh38,  int threads, string chemericOutJunction, string[] fastq_files, string outdir)
        {
            if (!GRCh37 && !GRCh38) return;

            if (GRCh37 && !Directory.Exists(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")))
                install(bin_directory, true, false);
            if (GRCh38 && !Directory.Exists(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play")))
                install(bin_directory, false, true);

            string read_command = fastq_files.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastq_files.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --annotate --examine_coding_effects" +
                (fastq_files.Length > 0 ? " --left_fq " + fastq_files[0] : "") +
                (fastq_files.Length > 1 ? " --right_fq " + fastq_files[1] : "") +
                " --CPU " + threads.ToString() +
                " --output_dir " + WrapperUtility.convert_windows_path(outdir) +
                " --genome_lib_dir " +
                    (GRCh37 ? WrapperUtility.convert_windows_path(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) :
                        WrapperUtility.convert_windows_path(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"))) +
                " --chimeric_junction " + WrapperUtility.convert_windows_path(chemericOutJunction);

            string script_name = Path.Combine(bin_directory, "star_fusion.bash");
            WrapperUtility.generate_and_run_script(script_name, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(bin_directory),
                "STAR/source/STAR-Fusion" + arguments
            }).WaitForExit();
            File.Delete(script_name);
        }

        public static void install(string current_directory, bool GRCh37, bool GRCh38)
        {
            string script_path = Path.Combine(current_directory, "install_star.bash");
            WrapperUtility.generate_and_run_script(script_path, new List<string>
            {
                "cd " + WrapperUtility.convert_windows_path(current_directory),
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "git clone https://github.com/alexdobin/STAR.git" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "cd STAR/source" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "make STAR" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "cd .." : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "export PATH=$PATH:" + WrapperUtility.convert_windows_path(Path.Combine(current_directory, "STAR", "source")) : "",
                //"perl -MCPAN -e 'install DB_File'",
                //"perl -MCPAN -e 'install URI::Escape'",
                //"perl -MCPAN -e 'install Set::IntervalTree'",
                //"perl -MCPAN -e 'install Carp::Assert'",
                //"perl -MCPAN -e 'install JSON::XS'",
                !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0")) ? "wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/v1.1.0/STAR-Fusion_v1.1.0.tar.gz" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0")) ? "tar -xvf STAR-Fusion_v1.1.0.tar.gz" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0")) ? "rm STAR-Fusion_v1.1.0.tar.gz" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0")) ? "cd STAR-Fusion_v1.1.0" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0")) ? "make" : "",
                "cd STAR-Fusion_v1.1.0",
                GRCh37 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) 
                ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                GRCh37 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) 
                ? "tar -xvf GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                GRCh37 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) 
                ? "rm GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                GRCh38 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play")) 
                ? "wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                GRCh38 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play")) 
                ? "tar -xvf GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
                GRCh38 && !Directory.Exists(Path.Combine(current_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play")) 
                ? "rm GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play.tar.gz" : "",
            }).WaitForExit();
            File.Delete(script_path);
        }
    }
}
