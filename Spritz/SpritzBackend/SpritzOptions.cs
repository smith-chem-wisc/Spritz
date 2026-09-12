using System;
using System.IO;
using System.Reflection;

namespace SpritzBackend
{
    public class SpritzOptions
    {
        public string AnalysisDirectory { get; set; }
        public string Fastq1 { get; set; }
        public string Fastq2 { get; set; }
        public string Fastq1SingleEnd { get; set; }
        public string SraAccession { get; set; }
        public string SraAccessionSingleEnd { get; set; }

        /// <summary>
        /// Filename of a user-supplied VCF in the analysis directory, used instead of calling variants
        /// from reads. A filename, not a path, for the same reason as Fastq1: only the analysis and
        /// resources directories are bind-mounted into the container.
        /// </summary>
        public string Vcf { get; set; }

        /// <summary>
        /// "vertebrates" (ftp.ensembl.org) or "bacteria" (Ensembl Genomes). Bacteria are a separate
        /// site with separate release numbering, so this is not derivable from the species name.
        /// </summary>
        public string Division { get; set; }

        /// <summary>
        /// "auto" (default), "ensembl" or "bootstrap" - where GATK base recalibration gets its
        /// known variant sites. "auto" asks Ensembl whether it publishes any for this species.
        /// </summary>
        public string KnownSites { get; set; }

        public int Threads { get; set; }

        /// <summary>podman (default), docker, or apptainer. See ContainerRuntime.</summary>
        public string ContainerRuntime { get; set; }
        public string Reference { get; set; }
        public bool AnalyzeVariants { get; set; }
        public bool AnalyzeIsoforms { get; set; }
        public bool Quantify { get; set; }
        public bool AvailableReferences { get; set; }
        public bool FetchGenomes { get; set; }
        public bool AnalysisSetup { get; set; }

        public static string DefaultAnalysisDirectory()
        {
            string defaultDirectory = Path.Combine(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location), "results");
            if (!RunnerEngine.IsDirectoryWritable(defaultDirectory))
            {
                defaultDirectory = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile), "Spritz", "output");
            }
            return defaultDirectory;
        }
    }
}