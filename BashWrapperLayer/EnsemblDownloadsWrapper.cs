using System;
using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    public class EnsemblDownloadsWrapper
    {

        private static string GRCh37PrimaryAssemblyUrl = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz";
        public static string GRCh37PrimaryAssemblyFilename { get; } = "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
        private static string GRCh38PrimaryAssemblyUrl = "ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
        public static string GRCh38PrimaryAssemblyFilename { get; } = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";

        private static string GRCh37GtfGeneModelUrl = "ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz";
        public static string GRCh37GtfGeneModelFilename { get; } = "Homo_sapiens.GRCh37.75.gtf";
        private static string GRCh38GtfGeneModelUrl = "ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz";
        public static string GRCh38GtfGeneModelFilename { get; } = "Homo_sapiens.GRCh38.81.gtf";

        // There are no gff3 gene models for GRCh37
        private static string GRCh38Gff3GeneModelUrl = "ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz";
        public static string GRCh38Gff3GeneModelFilename { get; } = "Homo_sapiens.GRCh38.81.gff3";

        public static void DownloadReferences(string binDirectory, string targetDirectory, string reference, out string genomeFastaPath, out string gtfGeneModelPath, out string gff3GeneModelPath)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);

            genomeFastaPath = downloadGrch37 ?
                Path.Combine(targetDirectory, GRCh37PrimaryAssemblyFilename) :
                downloadGrch38 ?
                    Path.Combine(targetDirectory, GRCh38PrimaryAssemblyFilename) :
                    "";
            gtfGeneModelPath = downloadGrch37 ?
                Path.Combine(targetDirectory, GRCh37GtfGeneModelFilename) :
                downloadGrch38 ?
                    Path.Combine(targetDirectory, GRCh38GtfGeneModelFilename) :
                    "";
            gff3GeneModelPath = downloadGrch37 ?
                gtfGeneModelPath :
                downloadGrch38 ?
                    Path.Combine(targetDirectory, GRCh38Gff3GeneModelFilename) :
                    "";

            if (!downloadGrch37 && !downloadGrch38)
                return;

            string scriptPath = Path.Combine(binDirectory, "scripts", "downloadEnsemblReference.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(targetDirectory),
                "if [ ! -f " + Path.GetFileName(genomeFastaPath) + " ]; then wget " + (downloadGrch38 ? GRCh38PrimaryAssemblyUrl : GRCh37PrimaryAssemblyUrl) + "; fi",
                "if [ -f " + Path.GetFileName(genomeFastaPath) + ".gz ]; then gunzip " + Path.GetFileName(genomeFastaPath) + ".gz; fi",
                "if [ ! -f " + Path.GetFileName(gtfGeneModelPath) + " ]; then wget " + (downloadGrch38 ? GRCh38GtfGeneModelUrl : GRCh37GtfGeneModelUrl) + "; fi",
                "if [ -f " + Path.GetFileName(gtfGeneModelPath) + ".gz ]; then gunzip " + Path.GetFileName(gtfGeneModelPath) + ".gz; fi",
                "if [ ! -f " + Path.GetFileName(gff3GeneModelPath) + " ]; then wget " + (downloadGrch38 ? GRCh38Gff3GeneModelUrl : GRCh37GtfGeneModelUrl) + "; fi", // note GRCh37 calls the gtf url instead
                "if [ -f " + Path.GetFileName(gff3GeneModelPath) + ".gz ]; then gunzip " + Path.GetFileName(gff3GeneModelPath) + ".gz; fi",
            }).WaitForExit();
        }
    }
}
