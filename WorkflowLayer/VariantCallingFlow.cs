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
        public List<string> GatkAnnotatedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> GatlSnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> GatkAnnotatedGenesSummaryPaths { get; private set; } = new List<string>();
        public List<string> ScapelVcfFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelAnnotatedVcfFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelSnpEffHtmlFilePaths { get; private set; } = new List<string>();
        public List<string> ScalpelAnnotatedGenesSummaryPaths { get; private set; } = new List<string>();

        public void CallVariants(string spritzDirectory, string reference, int threads, string genomeFasta, string sortedBed12Path, string ensemblKnownSitesPath,
            List<string> dedupedBamFiles, string reorderedFastaPath)
        {
            List<string> variantCallingCommands = new List<string>();
            string scriptName = Path.Combine(spritzDirectory, "scripts", "variantCalling.bash");
            foreach (string dedupedBam in dedupedBamFiles)
            {
                // GATK
                variantCallingCommands.AddRange(GATKWrapper.SplitNCigarReads(spritzDirectory, genomeFasta, dedupedBam, out string splitTrimBam));
                variantCallingCommands.AddRange(GATKWrapper.VariantCalling(spritzDirectory, threads, reorderedFastaPath, splitTrimBam, Path.Combine(spritzDirectory, ensemblKnownSitesPath), out string vcfPath));
                GatkVcfFilePaths.Add(vcfPath);
                variantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(spritzDirectory, reference, vcfPath, out string htmlReport, out string annotatedVcfPath, out string annotatedGenesSummaryPath));
                GatkAnnotatedVcfFilePaths.Add(annotatedVcfPath);
                GatlSnpEffHtmlFilePaths.Add(htmlReport);
                GatkAnnotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath);

                // Scalpel
                variantCallingCommands.AddRange(ScalpelWrapper.CallIndels(spritzDirectory, threads, genomeFasta, sortedBed12Path, splitTrimBam, Path.Combine(Path.GetDirectoryName(splitTrimBam), Path.GetFileNameWithoutExtension(splitTrimBam) + "_scalpelOut"), out string scalpelVcf));
                ScapelVcfFilePaths.Add(scalpelVcf);
                variantCallingCommands.AddRange(SnpEffWrapper.PrimaryVariantAnnotation(spritzDirectory, reference, scalpelVcf, out string htmlReport2, out string annotatedScalpelVcfPath, out string annotatedGenesSummaryPath2));
                ScalpelAnnotatedVcfFilePaths.Add(annotatedScalpelVcfPath);
                ScalpelSnpEffHtmlFilePaths.Add(htmlReport2);
                ScalpelAnnotatedGenesSummaryPaths.Add(annotatedGenesSummaryPath2);
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