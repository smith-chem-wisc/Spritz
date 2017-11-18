using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;

namespace RNASeqAnalysisWrappers
{
    public class STARWrapper
    {
        public static string BAM_SUFFIX = "Aligned.out.bam";

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
            string arguments = " --genomeLoad LoadAndExit" +
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
                "STAR/source/STAR --genomeLoad Remove"
            }).WaitForExit();
        }

        // fastqs must have \n line endings, not \r\n
        public static void BasicAlignReads(string bin_directory, int threads, string genomeDir, string[] fastq_files, string outprefix, bool strand_specific = true, string genomeLoad = "NoSharedMemory")
        {
            string reads_in = "\"" + String.Join("\" \"", fastq_files.Select(f => WrapperUtility.ConvertWindowsPath(f))) + "\"";
            string read_command = fastq_files.Any(f => Path.GetExtension(f) == ".gz") ?
                " --readFilesCommand zcat -c" :
                fastq_files.Any(f => Path.GetExtension(f) == ".bz2") ?
                    " --readFilesCommand bzip2 -c" :
                    "";
            string arguments =
                " --genomeLoad " + genomeLoad +
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
                "STAR/source/STAR" + arguments
            }).WaitForExit();
        }

        public void SuccessiveAlignReads(string bin_directory, int threads, string genomeDir, List<string[]> fastq_files, string outprefix, bool strand_specific = true)
        {
            LoadGenome(bin_directory, genomeDir);
            foreach(var fqs in fastq_files)
            {
                BasicAlignReads(bin_directory, threads, genomeDir, fqs, outprefix, strand_specific, "LoadAndKeep");
            }
            RemoveGenome(bin_directory);
        }

        public static void SubsetFastqs(string[] fastq_files, int reads, string current_directory, out string[] new_files)
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
        public static void RunStarFusion(string bin_directory, bool GRCh37, bool GRCh38,  int threads, string chemericOutJunction, string[] fastq_files, string outdir)
        {
            if (!GRCh37 && !GRCh38) return;

            if (GRCh37 && !Directory.Exists(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")))
                Install(bin_directory, true, false);
            if (GRCh38 && !Directory.Exists(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play")))
                Install(bin_directory, false, true);

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
                " --output_dir " + WrapperUtility.ConvertWindowsPath(outdir) +
                " --genome_lib_dir " +
                    (GRCh37 ? WrapperUtility.ConvertWindowsPath(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh37_gencode_v19_CTAT_lib_July192017.plug-n-play")) :
                        WrapperUtility.ConvertWindowsPath(Path.Combine(bin_directory, "STAR-Fusion_v1.1.0", "GRCh38_gencode_v26_CTAT_lib_July192017.plug-n-play"))) +
                " --chimeric_junction " + WrapperUtility.ConvertWindowsPath(chemericOutJunction);

            string script_name = Path.Combine(bin_directory, "scripts", "star_fusion.bash");
            WrapperUtility.GenerateAndRunScript(script_name, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(bin_directory),
                "STAR/source/STAR-Fusion" + arguments
            }).WaitForExit();
        }

        public static void Install(string current_directory, bool GRCh37, bool GRCh38)
        {
            string script_path = Path.Combine(current_directory, "install_star.bash");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(current_directory),
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "git clone https://github.com/alexdobin/STAR.git" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "cd STAR/source" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "make STAR" : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "cd .." : "",
                !Directory.Exists(Path.Combine(current_directory, "STAR")) ? "export PATH=$PATH:" + WrapperUtility.ConvertWindowsPath(Path.Combine(current_directory, "STAR", "source")) : "",
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
