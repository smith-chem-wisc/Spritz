using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;

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
        public int Workers { get; }
        public string DatabaseListPath { get; private set; }
        public string HtmlReportPath { get; private set; }
        public string AnnotatedVcfPath { get; private set; }
        public string AnnotatedGenesSummaryPath { get; private set; }
        public string VariantProteinFastaPath { get; private set; }
        public string VariantProteinXmlPath { get; private set; }

        public SnpEffWrapper(int workers)
        {
            Workers = workers;
        }

        /// <summary>
        /// Writes a script for installing SnpEff.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallSnpEff.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
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

        public List<string> PrimaryVariantAnnotation(string spritzDirectory, string reference, string inputVcfPath, bool fromReference = false)
        {
            string outPrefix = Path.Combine(Path.GetDirectoryName(inputVcfPath), Path.GetFileNameWithoutExtension(inputVcfPath));
            AnnotatedVcfPath = outPrefix + ".snpEffAnnotated.vcf";
            HtmlReportPath = outPrefix + ".snpEffAnnotated.html";
            AnnotatedGenesSummaryPath = outPrefix + ".snpEffAnnotated.genes.txt";
            VariantProteinFastaPath = outPrefix + ".snpEffAnnotated.protein.fasta";
            VariantProteinXmlPath = outPrefix + ".snpEffAnnotated.protein.xml";
            Directory.CreateDirectory(Path.Combine(spritzDirectory, "Tools", "SnpEff", "data"));
            string[] existingDatabases = Directory.GetDirectories(Path.Combine(spritzDirectory, "Tools", "SnpEff", "data"));
            if (File.Exists(AnnotatedVcfPath) && new FileInfo(AnnotatedVcfPath).Length > 0) return new List<string>();
            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                SnpEff(Workers) + " -v -stats " + WrapperUtility.ConvertWindowsPath(HtmlReportPath) +
                    " -fastaProt " + WrapperUtility.ConvertWindowsPath(VariantProteinFastaPath) +
                    " -xmlProt " + WrapperUtility.ConvertWindowsPath(VariantProteinXmlPath) +
                    " " + Path.GetFileName(existingDatabases.FirstOrDefault(x => Path.GetFileName(x).StartsWith(reference, true, null))) +
                    (fromReference ? "" : $" {WrapperUtility.ConvertWindowsPath(inputVcfPath)} > {WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath)}"),

                // ensure that the files get closed before continuing
                WrapperUtility.EnsureClosedFileCommands(WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath)),
                WrapperUtility.EnsureClosedFileCommands(WrapperUtility.ConvertWindowsPath(VariantProteinFastaPath)),
                WrapperUtility.EnsureClosedFileCommands(WrapperUtility.ConvertWindowsPath(VariantProteinXmlPath)),

                // remove the annotated VCF file if snpEff didn't work, e.g. if there was no VCF file to annotate
                "if [[ ( -f " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath) + " && ! -s " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath) + " ) ]]; then",
                "  rm " + WrapperUtility.ConvertWindowsPath(AnnotatedVcfPath),
                "fi",
            };
        }

        // see here for how to generate them from scratch: http://lab.loman.net/2012/11/16/how-to-get-snpeff-working-with-bacterial-genomes-from-ncbi/
        public void DownloadSnpEffDatabase(string spritzDirectory, string analysisDirectory, string reference)
        {
            DatabaseListPath = Path.Combine(spritzDirectory, "snpEffDatabases.txt");

            // check for existing list and database
            bool databaseListExists = File.Exists(DatabaseListPath);
            string databaseDirectory = Path.Combine(spritzDirectory, "Tools", "SnpEff", "data");
            string[] existingDatabases = Directory.Exists(databaseDirectory) ? Directory.GetDirectories(databaseDirectory) : new string[0];
            bool databaseExists = existingDatabases.Any(d => Path.GetFileName(d).StartsWith(reference, true, null));
            if (databaseListExists && databaseExists)
                return;

            // download database list
            string scriptPath = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "SnpEffDatabaseDownloadList.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "echo \"Downloading list of SnpEff references\"",
                SnpEff(Workers) + " databases > " + WrapperUtility.ConvertWindowsPath(DatabaseListPath),
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

            // download database (it downloads automatically now, with more feedback), but still need the mitochondrial references
            scriptPath = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "SnpEffDatabaseDownload.bash");
            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "echo \"\n# " + snpeffReference + "\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory,"Tools", "SnpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".genome : Human genome " + snpeffReference.Split('.')[0] + " using RefSeq transcripts\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory,"Tools", "SnpEff", "snpEff.config")),
                "echo \"" + snpeffReference + ".reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".M.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
                "echo \"\t" + snpeffReference + ".MT.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
            }).WaitForExit();
        }

        /// <summary>
        /// Gets string command for snpEff without arguments
        /// </summary>
        /// <returns></returns>
        public static string SnpEff(int workers)
        {
            var performance = new PerformanceCounter("Memory", "Available MBytes");
            var memory = performance.NextValue();
            return $"java -Xmx{(int)Math.Floor(memory / (double)workers)}M -jar SnpEff/snpEff.jar";
        }

        public static string GenerateXmlDatabaseFromReference(string spritzDirectory, string analysisDirectory, string reference, string inputFilePathForFilePrefix)
        {
            var snpeff = new SnpEffWrapper(1);
            WrapperUtility.GenerateAndRunScript(WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "SnpEffGenerateProteinXml.bash"), 
                snpeff.PrimaryVariantAnnotation(spritzDirectory, reference, inputFilePathForFilePrefix, true)).WaitForExit();
            return snpeff.VariantProteinXmlPath;
        }

        /// <summary>
        /// Creates a snpeff model for a custom gene model
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="genomeFastaPath"></param>
        /// <param name="geneModelGtfOrGffPath"></param>
        /// <returns>Name of the snpEff reference that was generated</returns>
        public static string GenerateDatabase(string spritzDirectory, string analysisDirectory, string genomeFastaPath, string referenceProteinFastaPath, string geneModelGtfOrGffPath)
        {
            string snpEffReferenceName = Path.GetExtension(geneModelGtfOrGffPath).Substring(1).ToUpperInvariant() + geneModelGtfOrGffPath.GetHashCode().ToString();
            string snpEffReferenceFolderPath = Path.Combine(spritzDirectory, "Tools", "SnpEff", "data", snpEffReferenceName);
            string scriptPath = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "SnpEffDatabaseGeneration.bash");
            string geneModelOption = Path.GetExtension(geneModelGtfOrGffPath).EndsWith("gtf") ? "-gtf22" : "-gff3";

            // if the database is already made, don't remake it
            if (File.Exists(Path.Combine(spritzDirectory, "Tools", "SnpEff", "data", snpEffReferenceName, "snpEffectPredictor.bin"))) { return snpEffReferenceName; }

            WrapperUtility.GenerateAndRunScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "cd SnpEff",

                // create data folder for this reference, and copy the custom gene model (can also copy regulatory annotations)
                "mkdir data/" + snpEffReferenceName,
                "cp " + WrapperUtility.ConvertWindowsPath(geneModelGtfOrGffPath) + " " + WrapperUtility.ConvertWindowsPath(Path.Combine(snpEffReferenceFolderPath, "genes" + Path.GetExtension(geneModelGtfOrGffPath))),
                "cp " + WrapperUtility.ConvertWindowsPath(referenceProteinFastaPath) + " " + WrapperUtility.ConvertWindowsPath(Path.Combine(snpEffReferenceFolderPath, "protein.fa")),

                // copy the genome to the genomes folder
                "mkdir " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "data", "genomes")),
                "cp " + WrapperUtility.ConvertWindowsPath(genomeFastaPath) + " " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "data", "genomes", snpEffReferenceName + ".fa")),

                // configure SnpEff for this custom reference
                // note: if different organism is used in the future, this becomes pretty complex... probably would list the organisms from snpEff.config in the GUI
                "echo \"\n# " + snpEffReferenceName + "\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
                "echo \"" + snpEffReferenceName + ".genome : Homo_sapiens\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory,"Tools", "SnpEff", "snpEff.config")),
                "echo \"" + snpEffReferenceName + ".reference : ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
                "echo \"\t" + snpEffReferenceName + ".M.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),
                "echo \"\t" + snpEffReferenceName + ".MT.codonTable : Vertebrate_Mitochondrial\" >> " + WrapperUtility.ConvertWindowsPath(Path.Combine(spritzDirectory, "Tools", "SnpEff", "snpEff.config")),

                // build snpEff model
                "cd ..",
                SnpEff(1) + " build " + geneModelOption + " -v " + snpEffReferenceName,
            }).WaitForExit();
            return snpEffReferenceName;
        }
    }
}