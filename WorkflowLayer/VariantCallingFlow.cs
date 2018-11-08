using Proteogenomics;
using System;
using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class VariantCallingFlow
    {
        /// <summary>
        /// These GVCF files contain phased information for every base in the genome. They're gigantic, so they're compressed by specifying .vcf.gz as the output.
        /// </summary>
        public List<string> GatkGvcfFilePaths { get; } = new List<string>();

        /// <summary>
        /// These are genotyped vcf files containing only variant information. They're the usable ones.
        /// </summary>
        public List<string> GatkVcfFilePaths { get; } = new List<string>();
        public List<string> GatkFilteredVcfFilePaths { get; } = new List<string>();
        public List<string> ScalpelVcfFilePaths { get; } = new List<string>();
        public List<string> ScalpelFilteredVcfFilePaths { get; } = new List<string>();
        public List<string> CombinedVcfFilePaths { get; } = new List<string>();
        public List<string> CombinedSortedVcfFilePaths { get; } = new List<string>();
        public List<string> CombinedAnnotatedVcfFilePaths { get; } = new List<string>();
        public List<string> CombinedSnpEffHtmlFilePaths { get; } = new List<string>();
        public List<string> CombinedAnnotatedGenesSummaryPaths { get; } = new List<string>();
        public List<string> CombinedAnnotatedProteinFastaPaths { get; } = new List<string>();
        public List<string> CombinedAnnotatedProteinXmlPaths { get; } = new List<string>();

        //public string CombinedGatkGvcfFilePath { get; private set; }
        //public string CombinedGatkVcfFilePath { get; private set; }
        //public string CombinedGatkFilteredVcfFilePath { get; private set; }
        //public string CombinedAnnotatedVcfFilePath { get; private set; }
        //public string CombinedSnpEffHtmlFilePath { get; private set; }
        //public string CombinedAnnotatedGenesSummaryPath { get; private set; }
        //public string CombinedAnnotatedProteinFastaPath { get; private set; }
        //public string CombinedAnnotatedProteinXmlPath { get; private set; }

        public void CallVariants(string spritzDirectory, string analysisDirectory, ExperimentType experimentType, string reference, int threads, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath, Genome genome, bool quickSnpEff, string indelFinder, int workers)
        {
            // Generate scripts for each BAM file
            List<string> variantCallingBashScripts = new List<string>();
            List<SnpEffWrapper> snpeffs = new List<SnpEffWrapper>();
            foreach (string dedupedBam in dedupedBamFiles)
            {
                List<string> variantCallingCommands = new List<string>();
                int workerThreads = (int)Math.Floor((double)threads / (double)workers);
                workerThreads = workerThreads == 0 ? workerThreads++ : workerThreads;

                // GATK
                var gatk = new GATKWrapper(workers);
                if (experimentType == ExperimentType.RNASequencing)
                {
                    variantCallingCommands.AddRange(gatk.SplitNCigarReads(spritzDirectory, reorderedFastaPath, dedupedBam));
                    variantCallingCommands.AddRange(gatk.BaseRecalibration(spritzDirectory, analysisDirectory, reorderedFastaPath, gatk.SplitTrimBamPath, ensemblKnownSitesPath));
                }
                else
                {
                    variantCallingCommands.AddRange(gatk.BaseRecalibration(spritzDirectory, analysisDirectory, reorderedFastaPath, dedupedBam, ensemblKnownSitesPath));
                }
                variantCallingCommands.AddRange(gatk.VariantCalling(spritzDirectory, experimentType, workerThreads, reorderedFastaPath, gatk.RecalibratedBamPath, Path.Combine(spritzDirectory, ensemblKnownSitesPath)));
                GatkGvcfFilePaths.Add(gatk.HaplotypeCallerGvcfPath);
                GatkVcfFilePaths.Add(gatk.HaplotypeCallerVcfPath);
                GatkFilteredVcfFilePaths.Add(gatk.FilteredHaplotypeCallerVcfPath);

                // Scalpel
                var scalpel = new ScalpelWrapper();
                bool useScalpel = indelFinder.Equals("scalpel", System.StringComparison.InvariantCultureIgnoreCase);
                if (useScalpel)
                {
                    variantCallingCommands.AddRange(scalpel.CallIndels(spritzDirectory, workerThreads, reorderedFastaPath, sortedBed12Path, dedupedBam, Path.Combine(Path.GetDirectoryName(dedupedBam), Path.GetFileNameWithoutExtension(dedupedBam) + "_scalpelOut")));
                    ScalpelVcfFilePaths.Add(scalpel.IndelVcfPath);
                    ScalpelFilteredVcfFilePaths.Add(scalpel.FilteredIndelVcfPath);
                }

                // Combine & Annotate
                var vcftools = new VcfToolsWrapper();
                var snpEff = new SnpEffWrapper(workers);
                var outprefix = Path.Combine(Path.GetDirectoryName(gatk.RecalibratedBamPath), Path.GetFileNameWithoutExtension(gatk.RecalibratedBamPath));
                if (useScalpel)
                {
                    variantCallingCommands.Add(vcftools.Concatenate(spritzDirectory, new string[] { gatk.FilteredHaplotypeCallerVcfPath, scalpel.FilteredIndelVcfPath }, outprefix));
                    variantCallingCommands.AddRange(gatk.SortVCF(spritzDirectory, analysisDirectory, vcftools.VcfConcatenatedPath, reorderedFastaPath));
                    CombinedVcfFilePaths.Add(vcftools.VcfConcatenatedPath);
                    CombinedSortedVcfFilePaths.Add(gatk.SortedVcfPath);
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.SortedVcfPath));
                }
                else if (indelFinder.Equals("gatk", System.StringComparison.InvariantCultureIgnoreCase))
                {
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.HaplotypeCallerVcfPath));
                }
                else
                {
                    variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.FilteredHaplotypeCallerVcfPath));
                }
                CombinedAnnotatedVcfFilePaths.Add(snpEff.AnnotatedVcfPath);
                CombinedSnpEffHtmlFilePaths.Add(snpEff.HtmlReportPath);
                CombinedAnnotatedProteinFastaPaths.Add(snpEff.VariantProteinFastaPath);
                CombinedAnnotatedProteinXmlPaths.Add(snpEff.VariantProteinXmlPath);
                snpeffs.Add(snpEff);

                string littleScriptName = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, $"VariantCalling{dedupedBam.GetHashCode().ToString()}.bash");
                WrapperUtility.GenerateScript(littleScriptName, variantCallingCommands);
                variantCallingBashScripts.Add(littleScriptName);
            }

            // Run the scripts in parallel
            string scriptName = WrapperUtility.GetAnalysisScriptPath(analysisDirectory, $"VariantCalling.bash");
            List<string> runnerCommands = new List<string>();
            List<int> runners = new List<int>();
            for (int i = 1; i <= variantCallingBashScripts.Count; i++)
            {
                // runs in parallel unless it's spawning enough workers or at the end of the line
                string logPath = Path.Combine(Path.GetDirectoryName(dedupedBamFiles[i - 1]), Path.GetFileNameWithoutExtension(dedupedBamFiles[i - 1]) + ".variantCalling.log");
                bool waitForWorkersToFinish = i % workers == 0 || i == variantCallingBashScripts.Count;
                runnerCommands.Add($"echo \"Running {variantCallingBashScripts[i - 1]} in the background. See {WrapperUtility.ConvertWindowsPath(logPath).Trim('"')} for output.\"");
                runnerCommands.Add($"bash {WrapperUtility.ConvertWindowsPath(variantCallingBashScripts[i - 1])} &> {WrapperUtility.ConvertWindowsPath(logPath)} &");
                runnerCommands.Add($"proc{i.ToString()}=$!");
                runners.Add(i);

                if (waitForWorkersToFinish)
                {
                    runners.ForEach(r => runnerCommands.Add($"wait $proc{r.ToString()}"));
                    runners.Clear();
                }
            }
            WrapperUtility.GenerateAndRunScript(scriptName, runnerCommands).WaitForExit();

            // Combine GVCFs and make a final database
            // This doesn't work because CombineGVCFs doesn't handle MNPs... MergeVcfs doesn't actually decide anything about overlapping variants... 
            // https://github.com/broadinstitute/gatk/issues/1385
            // Note: I'm not going to figure out how to merge scalpel VCFs just yet until I find out whether it does any better than GATK...
            //if (GatkVcfFilePaths.Count > 1 && !indelFinder.Equals("scalpel", System.StringComparison.InvariantCultureIgnoreCase))
            //{
            //    var gatk = new GATKWrapper();
            //    var snpEff = new SnpEffWrapper();
            //    variantCallingCommands.AddRange(gatk.CombineAndGenotypeGvcfs(spritzDirectory, reorderedFastaPath, GatkVcfFilePaths));
            //    if (indelFinder.Equals("gatk", System.StringComparison.InvariantCultureIgnoreCase))
            //    {
            //        variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.HaplotypeCallerVcfPath, quickSnpEff));
            //    }
            //    else
            //    {
            //        variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.FilteredHaplotypeCallerVcfPath, quickSnpEff));
            //    }
            //    CombinedGatkGvcfFilePath = gatk.HaplotypeCallerGvcfPath;
            //    CombinedGatkVcfFilePath = gatk.HaplotypeCallerVcfPath;
            //    CombinedAnnotatedVcfFilePath = snpEff.AnnotatedVcfPath;
            //    CombinedSnpEffHtmlFilePath = snpEff.HtmlReportPath;
            //    CombinedAnnotatedProteinFastaPath = snpEff.VariantProteinFastaPath;
            //    CombinedAnnotatedProteinXmlPath = snpEff.VariantProteinXmlPath;
            //}
        }
    }
}