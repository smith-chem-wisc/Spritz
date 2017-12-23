using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using Proteomics;

namespace ToolWrapperLayer
{
    public class EnsemblDownloadsWrapper
    {

        #region Primary Assembly Genome Fasta URLs and Filenames

        //ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/README
        private static string GRCh37PrimaryAssemblyUrl = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz";
        public static string GRCh37PrimaryAssemblyFilename { get; } = "Homo_sapiens.GRCh37.75.dna.primary_assembly.fa";
        //ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/README
        private static string GRCh38PrimaryAssemblyUrl = "ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz";
        public static string GRCh38PrimaryAssemblyFilename { get; } = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";

        #endregion Primary Assembly Genome Fasta URLs and Filenames

        #region GTF Gene Model URLs and Filenames

        private static string GRCh37GtfGeneModelUrl = "ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz";
        public static string GRCh37GtfGeneModelFilename { get; } = "Homo_sapiens.GRCh37.75.gtf";
        private static string GRCh38GtfGeneModelUrl = "ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz";
        public static string GRCh38GtfGeneModelFilename { get; } = "Homo_sapiens.GRCh38.81.gtf";

        #endregion GTF Gene Model URLs and Filenames

        #region GFF3 Gene Model URLs and Filenames

        // There are no gff3 gene models for GRCh37
        private static string GRCh38Gff3GeneModelUrl = "ftp://ftp.ensembl.org/pub/release-81/gff3/homo_sapiens/Homo_sapiens.GRCh38.81.gff3.gz";
        public static string GRCh38Gff3GeneModelFilename { get; } = "Homo_sapiens.GRCh38.81.gff3";

        #endregion GFF3 Gene Model URLs and Filenames

        #region Protein Fasta URLs and Filenames

        private static string GRCh37ProteinFastaUrl = "ftp://ftp.ensembl.org/pub/release-75//fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.75.pep.all.fa.gz";
        public static string GRCh37ProteinFastaFilename { get; } = "Homo_sapiens.GRCh37.75.pep.all.fa";
        private static string GRCh38ProteinFastaUrl = "ftp://ftp.ensembl.org/pub/release-81//fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz";
        public static string GRCh38ProteinFastaFilename { get; } = "Homo_sapiens.GRCh38.pep.all.fa";

        #endregion Protein Fasta URLs and Filenames

        public static void DownloadReferences(string binDirectory, string targetDirectory, string reference, out string genomeFastaPath, out string gtfGeneModelPath, out string gff3GeneModelPath, out string proteinFastaPath)
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
            proteinFastaPath = downloadGrch37 ?
                Path.Combine(targetDirectory, GRCh37ProteinFastaFilename) :
                downloadGrch38 ?
                    Path.Combine(targetDirectory, GRCh38ProteinFastaFilename) :
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
                "if [ ! -f " + Path.GetFileName(proteinFastaPath) + " ]; then wget " + (downloadGrch38 ? GRCh38ProteinFastaUrl : GRCh37ProteinFastaUrl) + "; fi", // note GRCh37 calls the gtf url instead
                "if [ -f " + Path.GetFileName(proteinFastaPath) + ".gz ]; then gunzip " + Path.GetFileName(proteinFastaPath) + ".gz; fi",
            }).WaitForExit();

            //Genome.WriteFasta(new Genome(genomeFastaPath).KaryotypicOrder(), genomeFastaPath); // todo: try this for ordering contigs before alignments; does gtf then need to be reordered?
        }

        public static Dictionary<string, string> UCSC2EnsemblChromosomeMappings(string binDirectory, string reference)
        {
            bool useGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool useGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            Dictionary<string, string> chromMappings = File.ReadAllLines(useGrch37 ?
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh37_UCSC2ensembl.txt") :
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh38_UCSC2ensembl.txt"))
                .Select(line => line.Split('\t'))
                .Where(x => x.Length > 1)
                .ToDictionary(line => line[0], line => line[1]);
            return chromMappings;
        }

        public static Dictionary<string, string> Ensembl2UCSCChromosomeMappings(string binDirectory, string reference)
        {
            bool useGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool useGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            Dictionary<string, string> chromMappings = File.ReadAllLines(useGrch37 ?
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh37_ensembl2UCSC.txt") :
                Path.Combine(binDirectory, "ChromosomeMappings", "GRCh38_ensembl2UCSC.txt"))
                .Select(line => line.Split('\t'))
                .Where(x => x.Length > 1)
                .ToDictionary(line => line[0], line => line[1]);
            return chromMappings;
        }

        /// <summary>
        /// Ensembl coding domain sequences (CDS) sometimes don't have start or stop codons annotated.
        /// The only way I can figure out how to tell which they are is to read in the protein FASTA and find the ones starting with X's or containing a stop codon '*'
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <param name="proteinFastaPath"></param>
        /// <returns></returns>
        public static void GetImportantProteinAccessions(string binDirectory, string proteinFastaPath, out HashSet<string> badProteinAccessions, out Dictionary<string, string> selenocysteineProteinAccessions)
        {
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(proteinFastaPath, true, DecoyType.None, false, ProteinDbLoader.ensembl_accession_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_fullName_expression, ProteinDbLoader.ensembl_gene_expression);
            HashSet<string> badOnes = new HashSet<string>(proteins.Where(p => p.BaseSequence.Contains('X') || p.BaseSequence.Contains('*')).Select(p => p.Accession));
            badProteinAccessions = badOnes;
            selenocysteineProteinAccessions = proteins.Where(p => !badOnes.Contains(p.Accession) && p.BaseSequence.Contains('U')).ToDictionary(p => p.Accession, p => p.BaseSequence);
        }
    }
}
