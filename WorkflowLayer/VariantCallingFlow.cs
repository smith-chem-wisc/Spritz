using Proteogenomics;
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
        public List<string> ScalpelFilteredVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedSortedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedSnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedGenesSummaryPaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedProteinFastaPaths { get; private set; } = new List<string>();
        public List<string> CombinedAnnotatedProteinXmlPaths { get; private set; } = new List<string>();

        public void CallVariants(string spritzDirectory, string analysisDirectory, string reference, int threads, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath, Genome genome, bool quickSnpEff, string indelFinder)
        {
            List<string> variantCallingCommands = new List<string>();
            List<SnpEffWrapper> snpeffs = new List<SnpEffWrapper>();
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
                bool useScalpel = indelFinder.Equals("scalpel", System.StringComparison.InvariantCultureIgnoreCase);
                if (useScalpel)
                {
                    variantCallingCommands.AddRange(scalpel.CallIndels(spritzDirectory, threads, reorderedFastaPath, sortedBed12Path, dedupedBam, Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam) + "_scalpelOut")));
                    ScalpelVcfFilePaths.Add(scalpel.IndelVcfPath);
                    ScalpelFilteredVcfFilePaths.Add(scalpel.FilteredIndelVcfPath);
                }

                // Combine & Annotate
                var vcftools = new VcfToolsWrapper();
                var snpEff = new SnpEffWrapper();
                var outprefix = Path.Combine(Path.GetDirectoryName(gatk.SplitTrimBamPath), Path.GetFileNameWithoutExtension(gatk.SplitTrimBamPath));
                if (useScalpel)
                {
                    variantCallingCommands.Add(vcftools.Concatenate(spritzDirectory, new string[] { gatk.FilteredHaplotypeCallerVcfPath, scalpel.FilteredIndelVcfPath }, outprefix));
                    variantCallingCommands.AddRange(gatk.SortVCF(spritzDirectory, analysisDirectory, vcftools.VcfConcatenatedPath, reorderedFastaPath));
                    CombinedVcfFilePaths.Add(vcftools.VcfConcatenatedPath);
                    CombinedSortedVcfFilePaths.Add(gatk.SortedVcfPath);
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, analysisDirectory, reference, gatk.SortedVcfPath, quickSnpEff));
                }
                else if (indelFinder.Equals("gatk", System.StringComparison.InvariantCultureIgnoreCase))
                {
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, analysisDirectory, reference, gatk.HaplotypeCallerVcfPath, quickSnpEff));
                }
                else
                {
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, analysisDirectory, reference, gatk.FilteredHaplotypeCallerVcfPath, quickSnpEff));
                }
                CombinedAnnotatedVcfFilePaths.Add(snpEff.AnnotatedVcfPath);
                CombinedSnpEffHtmlFilePaths.Add(snpEff.HtmlReportPath);
                CombinedAnnotatedProteinFastaPaths.Add(snpEff.VariantProteinFastaPath);
                CombinedAnnotatedProteinXmlPaths.Add(snpEff.VariantProteinXmlPath);
                snpeffs.Add(snpEff);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, variantCallingCommands).WaitForExit();
        }
    }
}