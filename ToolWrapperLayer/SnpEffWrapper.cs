using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

namespace ToolWrapperLayer
{
    /// <summary>
    /// SnpEff is a program for annotating variants, e.g. as missense or synonymous mutations.
    ///
    /// Program citation: http://snpeff.sourceforge.net/SnpEff_paper.pdf
    ///
    /// Note: SnpEff realigns indels to the 3-prime end, which isn't ideal based on the scalpel review.
    /// </summary>
    public class SnpEffWrapper :
        IInstallable
    {
        public string DatabaseListPath { get; private set; }
        public string HtmlReportPath { get; private set; }
        public string AnnotatedVcfPath { get; private set; }
        public string AnnotatedGenesSummaryPath { get; private set; }
        public string VariantProteinFastaPath { get; private set; }
        public string VariantProteinXmlPath { get; private set; }

        /// <summary>
        /// Writes a script for installing SnpEff.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "installScripts", "installSnpEff.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "if [ ! -d SnpEff ]; then git clone --depth=1 https://github.com/smith-chem-wisc/SnpEff; fi",
                "if [ ! -f SnpEff/snpEff.jar ]; then",
                "  cd SnpEff",
                "  mvn install:install-file -Dfile=lib/antlr-4.5.1-complete.jar -DgroupId=org.antlr -DartifactId=antlr -Dversion=4.5.1 -Dpackaging=jar",
                "  mvn install:install-file -Dfile=lib/biojava3-core-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-core -Dversion=3.0.7 -Dpackaging=jar",
                "  mvn install:install-file -Dfile=lib/biojava3-structure-3.0.7.jar -DgroupId=org.biojava -DartifactId=biojava3-structure -Dversion=3.0.7 -Dpackaging=jar",
                "  export VERSION=4.3",
                "  export VERSION_UND=`echo $VERSION | tr '.' '_'`",
                "  mvn clean compile assembly:assembly",
                "  mvn install:install-file -Dfile=target/SnpEff-$VERSION.jar -DgroupId=org.snpeff -DartifactId=SnpEff -Dversion=$VERSION -Dpackaging=jar -DgeneratePom=true --quiet",
                "  cp target/SnpEff-$VERSION-jar-with-dependencies.jar snpEff.jar",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing SnpEff.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        public List<string> PrimaryVariantAnnotation(string spritzDirectory, string reference, string vcfPath)
        {
            string outPrefix = Path.Combine(Path.GetDirectoryName(vcfPath), Path.GetFileNameWithoutExtension(vcfPath));
            AnnotatedVcfPath = outPrefix + ".snpEffAnnotated.vcf";
            HtmlReportPath = outPrefix + ".snpEffAnnotated.html";
            AnnotatedGenesSummaryPath = outPrefix + ".snpEffAnnotated.genes.txt";
            VariantProteinFastaPath = outPrefix + ".snpEffAnnotated.protein.fasta";
            VariantProteinXmlPath = outPrefix + ".snpEffAnnotated.protein.xml";
            string[] existingDatabases = Directory.GetDirectories(Path.Combine(spritzDirectory, "SnpEff", "data"));
            if (File.Exists(AnnotatedVcfPath)) return new List<string>();
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "snpEffAnnotation.bash");
            return new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                SnpEff() + " -v -stats " + WrapperUtility.ConvertWindowsPath(HtmlReportPath) +
                    " -fastaProt " + WrapperUtility.ConvertWindowsPath(VariantProteinFastaPath) +
                    " -xmlProt " + WrapperUtility.ConvertWindowsPath(VariantProteinXmlPath) +
                    " " + Path.GetFileName(existingDatabases.FirstOrDefault(x => Path.GetFileName(x).StartsWith(reference, true, null))) +
                    " " + WrapperUtility.ConvertWindowsPath(vcfPath) +
                    " > " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath),

                // remove the annotated VCF file if snpEff didn't work, e.g. if there was no VCF file to annotate
                "if [[ ( -f " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath) + " && ! -s " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath) + " ) ]]; then",
                "  rm " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath),
                "fi",
            };
        }

        public static string SnpEff()
        {
            var performance = new PerformanceCounter("Memory", "Available MBytes");
            var memory = performance.NextValue();
            return "java -Xmx" + Math.Floor(memory) + "M -jar SnpEff/snpEff.jar";
        }

        // see here for how to generate them from scratch: http://lab.loman.net/2012/11/16/how-to-get-snpeff-working-with-bacterial-genomes-from-ncbi/
        public void DownloadSnpEffDatabase(string spritzDirectory, string reference)
        {
            DatabaseListPath = Path.Combine(spritzDirectory, "snpEffDatabases.txt");

            // check for existing list and database
            bool databaseListExists = File.Exists(DatabaseListPath);
            string databaseDirectory = Path.Combine(spritzDirectory, "snpEff", "data");
            string[] existingDatabases = Directory.Exists(databaseDirectory) ? Directory.GetDirectories(databaseDirectory) : new string[0];
            bool databaseExists = existingDatabases.Any(d => Path.GetFileName(d).StartsWith(reference, true, null));
            if (databaseListExists && databaseExists)
                return;

            // download database list
            string scriptPath = Path.Combine(spritzDirectory, "scripts", "snpEffDatabaseDownload.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "echo \"Downloading list of SnpEff references\"",
                SnpEff() + " databases > " + WrapperUtility.ConvertWindowsPath(DatabaseListPath),
                WrapperUtility.EnsureClosedFileCommands(DatabaseListPath)
            }).WaitForExit();

            List<string> databases = new List<string>();
            using (StreamReader reader = new StreamReader(DatabaseListPath))
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
                "cd " + WrapperUtility.ConvertWindowsPath(spritzDirectory),
                "echo \"Downloading " + snpeffReference + " snpEff reference\"",
                SnpEff() + " download " + snpeffReference,
                "echo \"\n# " + snpeffReference + "\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "snpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".genome : Human genome " + snpeffReference.Split('.')[0] + " using RefSeq transcripts\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "snpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "snpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".M.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "snpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".MT.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "snpEff", "snpEff.config")),
            }).WaitForExit();
        }
    }
}