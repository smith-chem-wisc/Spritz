using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class VariantCallingFlow
    {
        public List<string> VariantCallingCommands { get; private set; } = new List<string>();
        public List<string> VcfFilePaths { get; private set; } = new List<string>();
        public List<string> AnnotatedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> SnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> AnnotatedGenesSummaryPaths { get; private set; } = new List<string>();
        public List<string> ScapelVcfFilePaths { get; private set; } = new List<string>();
        public List<string> AnnotatedScapelVcfFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelSnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelAnnotatedGenesSummaryPaths { get; private set; } = new List<string>();

        public void CallVariants(string binDirectory, string reference, int threads, string genomeFasta, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath)
        {
            string scriptName = Path.Combine(binDirectory, "scripts", "variantCalling.bash");

            foreach (string dedupedBam in dedupedBamFiles)
            {
                // GATK
                VariantCallingCommands.AddRange(GATKWrapper.SplitNCigarReads(binDirectory, genomeFasta, dedupedBam, out string splitTrimBam));
                VariantCallingCommands.AddRange(GATKWrapper.VariantCalling(binDirectory, threads, reorderedFastaPath, splitTrimBam, Path.Combine(binDirectory, ensemblKnownSitesPath), out string vcfPath));
                VcfFilePaths.Add(vcfPath);
                VariantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(binDirectory, reference, vcfPath, out string htmlReport, out string annotatedVcfPath, out string annotatedGenesSummaryPath));
                AnnotatedVcfFilePaths.Add(annotatedVcfPath);
                SnpEffHtmlFilePaths.Add(htmlReport);
                AnnotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath);

                // Scalpel
                VariantCallingCommands.AddRange(ScalpelWrapper.CallIndels(binDirectory, threads, genomeFasta, sortedBed12Path, splitTrimBam, Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + "_scalpelOut"), out string scalpelVcf));
                ScapelVcfFilePaths.Add(scalpelVcf);
                VariantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(binDirectory, reference, scalpelVcf, out string htmlReport2, out string annotatedScalpelVcfPath, out string annotatedGenesSummaryPath2));
                AnnotatedScapelVcfFilePaths.Add(annotatedScalpelVcfPath);
                ScalpelSnpEffHtmlFilePaths.Add(htmlReport2);
                ScalpelAnnotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath2);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, VariantCallingCommands).WaitForExit();
        }
    }
}