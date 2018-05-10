using System.Collections.Generic;
using System.IO;
using ToolWrapperLayer;
using Proteogenomics;
using Bio.VCF;
using System;
using System.Linq;

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

        public void CallVariants(string spritzDirectory, string reference, int threads, string genomeFasta, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath)
        {
            List<string> variantCallingCommands = new List<string>();
            string scriptName = Path.Combine(spritzDirectory, "scripts", "variantCalling.bash");
            foreach (string dedupedBam in dedupedBamFiles)
            {
                // GATK
                var gatk = new GATKWrapper();
                variantCallingCommands.AddRange(gatk.SplitNCigarReads(spritzDirectory, genomeFasta, dedupedBam));
                variantCallingCommands.AddRange(gatk.VariantCalling(spritzDirectory, threads, reorderedFastaPath, gatk.SplitTrimBamPath, Path.Combine(spritzDirectory, ensemblKnownSitesPath)));
                GatkVcfFilePaths.Add(gatk.HaplotypeCallerVcfPath);
                GatkFilteredVcfFilePaths.Add(gatk.FilteredHaplotypeCallerVcfPath);

                // Scalpel
                var scalpel = new ScalpelWrapper();
                variantCallingCommands.AddRange(scalpel.CallIndels(spritzDirectory, threads, genomeFasta, sortedBed12Path, gatk.SplitTrimBamPath, Path.Combine(Path.GetDirectoryName(gatk.SplitTrimBamPath), Path.GetFileNameWithoutExtension(gatk.SplitTrimBamPath) + "_scalpelOut")));
                ScalpelVcfFilePaths.Add(scalpel.IndelVcfPath);
                ScalpelFilteredlVcfFilePaths.Add(scalpel.FilteredIndelVcfPath);

                // Combine & Annotate
                var vcftools = new VcfToolsWrapper();
                var outprefix = Path.Combine(Path.GetDirectoryName(gatk.SplitTrimBamPath), Path.GetFileNameWithoutExtension(gatk.SplitTrimBamPath));
                variantCallingCommands.Add(vcftools.Concatenate(spritzDirectory, new string[] { gatk.FilteredHaplotypeCallerVcfPath, scalpel.FilteredIndelVcfPath }, outprefix));
                variantCallingCommands.AddRange(gatk.SortVCF(spritzDirectory, vcftools.VcfConcatenatedPath, genomeFasta));
                CombinedVcfFilePaths.Add(vcftools.VcfConcatenatedPath);
                CombinedSortedVcfFilePaths.Add(gatk.SortedVcfPath);
                var snpEff = new SnpEffWrapper();
                variantCallingCommands.AddRange(snpEff.PrimaryVariantAnnotation(spritzDirectory, reference, gatk.SortedVcfPath));
                CombinedAnnotatedVcfFilePaths.Add(snpEff.AnnotatedVcfPath);
                CombinedSnpEffHtmlFilePaths.Add(snpEff.HtmlReportPath);
                CombinedAnnotatedProteinFastaPaths.Add(snpEff.VariantProteinFastaPath);
                CombinedAnnotatedProteinXmlPaths.Add(snpEff.VariantProteinXmlPath);
            }
            WrapperUtility.GenerateAndRunScript(scriptName, variantCallingCommands).WaitForExit();
        }

        public List<Transcript> ApplyIndels(string indelVcfPath, Genome genome, GeneModel geneModel)
        {
            if (!File.Exists(indelVcfPath)) { throw new ArgumentException("Error: VCF not found: " + indelVcfPath); }

            // Parse Indels from VCF file
            VCFParser vcf = new VCFParser(indelVcfPath);
            List<Variant> indels = vcf.Select(x => x)
                .Where(x => x.AlternateAlleles.All(a => a.Length != x.Reference.Length))
                .Select(v => new Variant(null, v, genome)).ToList();
            return geneModel.ApplyVariants(indels);
        }

        // Code in and annotate 
        public List<Transcript> ApplySNVs(string snvVcfPath, Genome genome, GeneModel geneModel)
        {
            if (!File.Exists(snvVcfPath)) { throw new ArgumentException("Error: VCF not found: " + snvVcfPath); }

            // Parse SNVs from VCF file
            VCFParser vcf = new VCFParser(snvVcfPath);
            List<Variant> singleNucleotideVariants = vcf.Select(x => x)
                .Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1))
                .Select(v => new Variant(null, v, genome)).ToList();
            return geneModel.ApplyVariants(singleNucleotideVariants);
        }

        public void AnnotateSAVs(string snvVcfPath, Genome genome, GeneModel geneModel)
        {
            if (!File.Exists(snvVcfPath)) { throw new ArgumentException("Error: VCF not found: " + snvVcfPath); }

            // Parse SNVs from VCF file
            VCFParser vcf = new VCFParser(snvVcfPath);
            List<Variant> singleNucleotideVariants = vcf.Select(x => x)
                .Where(x => x.AlternateAlleles.All(a => a.Length == x.Reference.Length && a.Length == 1))
                .Select(v => new Variant(null, v, genome)).ToList();
            geneModel.AddVariantAnnotations(singleNucleotideVariants);
        }
    }
}