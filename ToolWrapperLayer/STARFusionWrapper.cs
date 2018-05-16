using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ToolWrapperLayer
{
    /// <summary>
    /// STARFusion is a program for analyzing gene fusion events detected using the STAR aligner.
    /// </summary>
    public class STARFusionWrapper :
        IInstallable
    {
        private readonly string Grch37SourceDataUrl = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_v19_CTAT_lib_Feb092018.source_data.tar.gz";
        private readonly string Grch37ReferenceUrl = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh37_v19_CTAT_lib_Feb092018.plug-n-play.tar.gz";
        private readonly string Grch38SourceDataUrl = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.source_data.tar.gz";
        private readonly string Grch38ReferenceUrl = "https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz";
        private readonly string PfamDataUrl = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz";

        public string ReferenceLibraryDirectory { get; private set; }
        public string OutputDirectoryPath { get; private set; }
        public string CodingEffectFilename { get; } = "star-fusion.fusion_predictions.abridged.coding_effect.tsv";

        /// <summary>
        /// Writes a script for installing STAR-Fusion.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallStarFusion.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d STAR-Fusion-v1.4.0 ]; then",
                "  wget https://github.com/STAR-Fusion/STAR-Fusion/releases/download/STAR-Fusion-v1.4.0/STAR-Fusion-v1.4.0.FULL.tar.gz",
                "  tar -xvf STAR-Fusion-v1.4.0.FULL.tar.gz",
                "  rm STAR-Fusion-v1.4.0.FULL.tar.gz",
                "  cd STAR-Fusion-v1.4.0",
                "  make",
                "  cd ..",

                // need hmmer for building STAR-Fusion references
                "  " + WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "  wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz",
                "  tar zxvf hmmer-3.1b2-linux-intel-x86_64.tar.gz",
                "  rm hmmer-3.1b2-linux-intel-x86_64.tar.gz",
                "  cd hmmer-3.1b2-linux-intel-x86_64",
                "  cp binaries/* /usr/local/bin",

                // jellyfish
                "  " + WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "  wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz",
                "  tar xvf jellyfish-2.2.7.tar.gz; rm jellyfish-2.2.7.tar.gz",
                "  cd jellyfish-2.2.7",
                "  ./configure; make; make install",

                // salmon
                "  " + WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "  wget https://github.com/COMBINE-lab/salmon/releases/download/v0.9.1/Salmon-0.9.1_linux_x86_64.tar.gz",
                "  tar xvf Salmon-0.9.1_linux_x86_64.tar.gz; rm Salmon-0.9.1_linux_x86_64.tar.gz",
                "  cp tar xvf Salmon-0.9.1_linux_x86_64/bin/salmon /usr/local/bin",
                "fi",
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing STAR-Fusion.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "RemoveStarFusion.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "rm -rf STAR-Fusion-v1.4.0"
            });
            return scriptPath;
        }

        /// <summary>
        /// Downloads precompiled reference (takes ~90 mins)
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="reference"></param>
        /// <returns></returns>
        public List<string> DownloadPrecompiledReference(string spritzDirectory, string reference)
        {
            bool grch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase);
            bool grch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase);
            string referenceUrl = grch37 ? Grch37ReferenceUrl : Grch38ReferenceUrl;
            string referenceTarGz = referenceUrl.Split('/').Last();
            string referenceFolderName = referenceTarGz.Split('.').First();
            ReferenceLibraryDirectory = Path.Combine(spritzDirectory, "Tools", "STAR-Fusion_v1.4.0", "data", referenceFolderName);

            bool downloadGrch37 = grch37 && !Directory.Exists(Path.Combine(spritzDirectory, "Tools", "STAR-Fusion_v1.4.0", "data", ReferenceLibraryDirectory));
            bool downloadGrch38 = grch38 && !Directory.Exists(Path.Combine(spritzDirectory, "Tools", "STAR-Fusion_v1.4.0", "data", ReferenceLibraryDirectory));

            if (!downloadGrch37 && !downloadGrch38) { return new List<string>(); }

            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "cd STAR-Fusion_v1.4.0",
                "mkdir data; cd data",
                "wget " + referenceUrl,
                "tar xvf " + referenceTarGz + "; rm " + referenceTarGz,
            };
        }

        /// <summary>
        /// Prepares reference from scratch (takes > 4 hours, and then I stopped it)
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="reference"></param>
        /// <param name="genomeFastaPath"></param>
        /// <param name="geneModelGtfPath"></param>
        /// <returns></returns>
        public List<string> PrepareReference(string spritzDirectory, int threads, string reference, string genomeFastaPath, string geneModelGtfPath)
        {
            bool downloadGrch37 = String.Equals(reference, "GRCh37", StringComparison.CurrentCultureIgnoreCase) &&
                !Directory.Exists(Path.Combine(spritzDirectory, "STAR-Fusion_v1.4.0", "data", "GRCh37_v19_CTAT_lib_Feb092018.source_data"));
            bool downloadGrch38 = String.Equals(reference, "GRCh38", StringComparison.CurrentCultureIgnoreCase) &&
                !Directory.Exists(Path.Combine(spritzDirectory, "STAR-Fusion_v1.4.0", "data", "GRCh38_v27_CTAT_lib_Feb092018.source_data"));

            if (!downloadGrch37 && !downloadGrch38) { return new List<string>(); }

            // CTAT resource library and reference genome
            string sourceDataUrl = downloadGrch37 ? Grch37SourceDataUrl : Grch38SourceDataUrl;
            string sourceDataTarGz = sourceDataUrl.Split('/').Last();
            string sourceDataFolderName = sourceDataTarGz.Split('.').First();
            ReferenceLibraryDirectory = Path.Combine(spritzDirectory, "Tools", "STAR-Fusion_v1.4.0", "data", sourceDataFolderName);

            // pfam data
            string pfamFilenameGz = PfamDataUrl.Split('/').Last();
            string pfamFilename = Path.GetFileNameWithoutExtension(pfamFilenameGz);

            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "cd STAR-Fusion_v1.4.0",
                "mkdir data; cd data",
                "wget " + PfamDataUrl + " " + sourceDataUrl,
                "gunzip Pfam-A.hmm.gz",
                "tar -xvf " + sourceDataTarGz + "; rm " + sourceDataTarGz,
                "gunzip " + pfamFilenameGz,
                "hmmpress " + pfamFilename,
                "cd ..",
                "FusionFilter/prep_genome_lib.pl " +
                    "--genome_fa " + WrapperUtility.ConvertWindowsPath(Path.Combine("data", sourceDataFolderName, "ref_genome.fa")) +
                    "--gtf " + WrapperUtility.ConvertWindowsPath(Path.Combine("data", sourceDataFolderName, "ref_annot.gtf")) +
                    "--fusion_annot_lib " + WrapperUtility.ConvertWindowsPath(Path.Combine("data", sourceDataFolderName, "CTAT_HumanFusionLib.v0.1.0.dat.gz")) +
                    "--annot_filter_rule " + WrapperUtility.ConvertWindowsPath(Path.Combine("data", sourceDataFolderName, "AnnotFilterRule.pm")) +
                    "--pfam_db " + WrapperUtility.ConvertWindowsPath(Path.Combine("data", pfamFilename)) +
                    "--CPU " + threads.ToString()
            };
        }

        /// <summary>
        /// Run star fusion from chimericOutJunctions
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="chimericOutJunction"></param>
        /// <param name="outdir"></param>
        /// <returns></returns>
        public List<string> RunStarFusion(string spritzDirectory, string analysisDirectory, int threads, string chimericOutJunction)
        {
            if (ReferenceLibraryDirectory == null)
            {
                throw new FileNotFoundException("STAR-Fusion reference library was not generated prior to running STAR-Fusion.");
            }

            OutputDirectoryPath = Path.Combine(Path.GetDirectoryName(chimericOutJunction), Path.GetFileNameWithoutExtension(chimericOutJunction) + "FusionAnalysis");
            Directory.CreateDirectory(OutputDirectoryPath);
            string tmp = Path.Combine(analysisDirectory, "_STARFusionTmp");
            Directory.CreateDirectory(tmp);

            string arguments =
                " --examine_coding_effects" +
                " --CPU " + threads.ToString() +
                " --output_dir " + WrapperUtility.ConvertWindowsPath(OutputDirectoryPath) +
                " --genome_lib_dir " + WrapperUtility.ConvertWindowsPath(ReferenceLibraryDirectory) +
                " --chimeric_junction " + WrapperUtility.ConvertWindowsPath(chimericOutJunction) +
                " --tmpdir " + WrapperUtility.ConvertWindowsPath(tmp);

            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "STAR-Fusion_v1.4.0/STAR-Fusion " + arguments
            };
        }

        /// <summary>
        /// Run star fusion from fastqs
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <param name="analysisDirectory"></param>
        /// <param name="threads"></param>
        /// <param name="fastqs"></param>
        /// <param name="outdir"></param>
        /// <returns></returns>
        public List<string> RunStarFusion(string spritzDirectory, string analysisDirectory, int threads, string[] fastqs)
        {
            if (ReferenceLibraryDirectory == null)
            {
                throw new FileNotFoundException("STAR-Fusion reference library was not generated prior to running STAR-Fusion.");
            }
            if (fastqs == null || fastqs.Length == 0)
            {
                throw new ArgumentException("No fastqs were passed into STAR-Fusion.");
            }

            OutputDirectoryPath = Path.Combine(Path.GetDirectoryName(fastqs[0]), Path.GetFileNameWithoutExtension(fastqs[0]) + "FusionAnalysis");
            Directory.CreateDirectory(OutputDirectoryPath);
            string tmp = Path.Combine(analysisDirectory, "_STARFusionTmp");
            Directory.CreateDirectory(tmp);

            string arguments =
                " --examine_coding_effects" +
                " --left_fq " + fastqs[0] +
                (fastqs.Length > 1 ? " --right_fq " + fastqs[1] : "") +
                " --CPU " + threads.ToString() +
                " --output_dir " + WrapperUtility.ConvertWindowsPath(OutputDirectoryPath) +
                " --genome_lib_dir " + WrapperUtility.ConvertWindowsPath(ReferenceLibraryDirectory) +
                " --tmpdir " + WrapperUtility.ConvertWindowsPath(tmp);

            return new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "STAR-Fusion_v1.4.0/STAR-Fusion" + arguments
            };
        }

        /// <summary>
        /// Takes a
        /// </summary>
        /// <param name="codingEffectFilePath"></param>
        /// <returns></returns>
        public List<Protein> ParseCodingEffect(string codingEffectFilePath, int minPeptideLength, string organism, HashSet<string> usedAccessions)
        {
            List<Protein> fusionProteins = new List<Protein>();
            using (StreamReader reader = new StreamReader(codingEffectFilePath))
            {
                while (true)
                {
                    string line = reader.ReadLine();
                    if (line == null) { break; }
                    if (line.StartsWith("#")) { break; } // skip header
                    string[] s = line.Split('\t');
                    var fusionName = s[0];
                    var junctionSpanningReadCount = s[1];
                    var spliceType = s[3];
                    var leftGene = s[4];
                    var leftBreakpoint = s[5];
                    var rightGene = s[6];
                    var rightBreakpoint = s[7];
                    var ffpm = s[9];
                    var annots = s[14];
                    var proteinFusionType = s[19];
                    var proteinFusionTranslation = s[22];
                    var proteinSequence = proteinFusionTranslation.Split('*')[0];
                    if (proteinSequence.Length < minPeptideLength) { continue; }
                    string accession = fusionName;
                    int i = 1;
                    while (usedAccessions.Contains(accession))
                    {
                        accession = fusionName + "_" + i++.ToString();
                    }
                    usedAccessions.Add(accession);
                    fusionProteins.Add(new Protein(proteinSequence, accession, organism, new List<Tuple<string, string>> { new Tuple<string, string>("fusion", leftGene + "--" + rightGene) }));
                }
            }
            return fusionProteins;
        }
    }
}