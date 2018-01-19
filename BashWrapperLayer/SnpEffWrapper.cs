using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Globalization;

namespace ToolWrapperLayer
{
    /// <summary>
    /// Program citation: http://snpeff.sourceforge.net/SnpEff_paper.pdf
    /// 
    /// Note: SnpEff realigns indels to the 3-prime end, which isn't ideal based on the scalpel review.
    /// </summary>
    public class SnpEffWrapper
    {

        #region Public Methods

        public static string SnpEff()
        {
            var performance = new PerformanceCounter("Memory", "Available MBytes");
            var memory = performance.NextValue();
            return "java -Xmx" + Math.Floor(memory) + "M -jar snpEff/snpEff.jar";
        }

        public static void PrimaryVariantAnnotation(string binDirectory, string reference, string vcfPath, out string htmlReport, out string annotatedVcfPath)
        {
            annotatedVcfPath = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath) + ".snpEffAnnotated.vcf");
            htmlReport = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath) + ".snpEffAnnotated.html");
            string[] existingDatabases = Directory.GetDirectories(Path.Combine(binDirectory, "snpEff", "data"));
            if (File.Exists(annotatedVcfPath)) return;
            string scriptPath = Path.Combine(binDirectory, "scripts", "snpEffAnnotation.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                SnpEff() + " -v -stats " + WrapperUtility.ConvertWindowsPath(htmlReport) +
                    " " + Path.GetFileName(existingDatabases.FirstOrDefault(x => Path.GetFileName(x).StartsWith(reference, true, null))) +
                    " " + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " > " + WrapperUtility.ConvertWindowsPath(annotatedVcfPath)
            }).WaitForExit();
        }

        // see here for how to generate them from scratch: http://lab.loman.net/2012/11/16/how-to-get-snpeff-working-with-bacterial-genomes-from-ncbi/
        public static void DownloadSnpEffDatabase(string binDirectory, string reference, out string databaseListPath)
        {
            databaseListPath = Path.Combine(binDirectory, "snpEffDatabases.txt");

            // check for existing list and database
            bool databaseListExists = File.Exists(databaseListPath);
            string databaseDirectory = Path.Combine(binDirectory, "snpEff", "data");
            string[] existingDatabases = Directory.Exists(databaseDirectory) ? Directory.GetDirectories(databaseDirectory) : new string[0];
            bool databaseExists = existingDatabases.Any(d => Path.GetFileName(d).StartsWith(reference, true, null));
            if (databaseListExists && databaseExists)
                return;

            // download database list
            string scriptPath = Path.Combine(binDirectory, "scripts", "snpEffDatabaseDownload.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "echo \"Downloading list of SnpEff references\"",
                SnpEff() + " databases > " + WrapperUtility.ConvertWindowsPath(databaseListPath),
                WrapperUtility.EnsureClosedFileCommands(databaseListPath)
            }).WaitForExit();

            List<string> databases = new List<string>();
            using (StreamReader reader = new StreamReader(databaseListPath))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null) break;
                    databases.Add(line.Split('\t')[0].TrimEnd());
                }
            }
            string snpeffReference = databases.FirstOrDefault(d => d.StartsWith(reference, true, CultureInfo.InvariantCulture));

            // download database
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "echo \"Downloading " + snpeffReference + " snpEff reference\"",
                SnpEff() + " download " + snpeffReference,
                "echo \"\n# " + snpeffReference + "\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "snpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".genome : Human genome " + snpeffReference.Split('.')[0] + " using RefSeq transcripts\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "snpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "snpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".M.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "snpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".MT.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(binDirectory, "snpEff", "snpEff.config")),
            }).WaitForExit();
        }

        public static void Install(string binDirectory)
        {
            if (Directory.Exists(Path.Combine(binDirectory, "snpEff")) && Directory.Exists(Path.Combine(binDirectory, "skewer-0.2.2")))
                return;
            string scriptPath = Path.Combine(binDirectory, "scripts", "snpEffInstaller.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(binDirectory),
                "wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip",
                "unzip snpEff_latest_core.zip",
                "rm snpEff_latest_core.zip",
            }).WaitForExit();
        }

        #endregion Public Methods

    }
}
