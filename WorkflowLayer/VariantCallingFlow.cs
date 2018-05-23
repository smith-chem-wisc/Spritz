using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class VariantCallingFlow
    {
        public List<string> GatkVcfFilePaths { get; private set; } = new List<string>();
        public List<string> GatkFilteredVcfFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelVcfFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelFilteredlVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedSortedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedSnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedGenesSummaryPaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedProteinFastaPaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedProteinXmlPaths { get; private set; } = new List<string>();

        public void CallVariants(string spritzDirectory, string analysisDirectory, string reference, int threads, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath)
        {
            new SnpEffWrapper().DownloadSnpEffDatabase(spritzDirectory, analysisDirectory, reference);
            List<string> variantCallingCommands = new List<string>();
            string scriptName = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, "VariantCalling.bash");
            foreach (string dedupedBam in dedupedBamFiles)
            {
                // GATK
                var gatk = new GATKWrapper();
                variantCallingCommands.AddRange(gatk.SplitNCigarReads(spritzDirectory, reorderedFastaPath, dedupedBam));
                variantCallingCommands.AddRange(gatk.VariantCalling(spritzDirectory, threads, reorderedFastaPath, gatk.SplitTrimBamPath, Path.Combine(spritzDirectory, ensemblKnownSitesPath)));
                GatkVcfFilePaths.Add(gatk.HaplotypeCallerVcfPath);
                GatkFilteredVcfFilePaths.Add(gatk.FilteredHaplotypeCallerVcfPath);

                // Scalpel
                var scalpel = new ScalpelWrapper();
                variantCallingCommands.AddRange(scalpel.CallIndels(spritzDirectory, threads, reorderedFastaPath, sortedBed12Path, dedupedBam, Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam) + "_scalpelOut")));
                ScalpelVcfFilePaths.Add(scalpel.IndelVcfPath);
                ScalpelFilteredlVcfFilePaths.Add(scalpel.FilteredIndelVcfPath);

                // Combine & Annotate
                var vcftools = new VcfToolsWrapper();
                var outprefix = Path.Combine(Path.GetDirectoryName(gatk.SplitTrimBamPath), Path.GetFileNameWithoutExtension(gatk.SplitTrimBamPath));
                variantCallingCommands.Add(vcftools.Concatenate(spritzDirectory, new string[] { gatk.FilteredHaplotypeCallerVcfPath, scalpel.FilteredIndelVcfPath }, outprefix));
                variantCallingCommands.AddRange(gatk.SortVCF(spritzDirectory, analysisDirectory, vcftools.VcfConcatenatedPath, reorderedFastaPath));
                CombinedVcfFilePaths.Add(vcftools.VcfConcatenatedPath);
                CombinedSortedVcfFilePaths.Add(gatk.SortedVcfPath);
                var snpEff = new SnpEffWrapper();
                variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, analysisDirectory, reference, gatk.SortedVcfPath));
                CombinedAnnotatedVcfFilePaths.Add(snpEff.AnnotatedVcfPath);
                CombinedSnpEffHtmlFilePaths.Add(snpEff.HtmlReportPath);
                CombinedAnnotatedProteinFastaPaths.Add(snpEff.VariantProteinFastaPath);
                CombinedAnnotatedProteinXmlPaths.Add(snpEff.VariantProteinXmlPath);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, variantCallingCommands).WaitForExit();
        }
    }
}