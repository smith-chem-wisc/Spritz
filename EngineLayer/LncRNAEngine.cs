using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using ToolWrapperLayer;
using UsefulProteomicsDatabases;
using System;

namespace WorkflowLayer
{
    public class LncRNAEngine
    {
        public static void RunLnckyRNAFromFastq(string bin, string analysisDirectory, string reference, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string proteinFasta, string geneModelGtfOrGff, string ensemblKnownSitesPath, bool useReadSubset = false, int readSubset = 300000)
        {
            if (Path.GetExtension(genomeFasta) == ".gz")
            {
                WrapperUtility.RunBashCommand("gunzip", WrapperUtility.ConvertWindowsPath(genomeFasta));
                genomeFasta = Path.ChangeExtension(genomeFasta, null);
            }

            // We need to use the same fasta file throughout and have all the VCF and GTF chromosome reference IDs be the same as these.
            // Right now this is based on ensembl references, so those are the chromosome IDs I will be using throughout
            // TODO: try this with UCSC references to judge whether there's a difference in quality / yield / FDR etc in subsequent proteomics analysis
            // This file needs to be in karyotypic order; this allows us not to have to reorder it for GATK analysis
            string ensemblFastaHeaderDelimeter = " ";
            string reorderedFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".karyotypic.fa");
            Genome ensemblGenome = new Genome(genomeFasta);
            if (!ensemblGenome.IsKaryotypic(ensemblFastaHeaderDelimeter))
            {
                ensemblGenome.Chromosomes = ensemblGenome.KaryotypicOrder(ensemblFastaHeaderDelimeter);
                if (!File.Exists(reorderedFasta))
                {
                    Genome.WriteFasta(ensemblGenome.Chromosomes, reorderedFasta);
                }
            }
            else
            {
                reorderedFasta = genomeFasta;
            }

            // Alignment preparation
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "genomeGenerate.bash"), STARWrapper.GenerateGenomeIndex(bin, threads, genomeStarIndexDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff)).WaitForExit();

            List<string[]> fastqsForAlignment = new List<string[]>();
            List<bool> strandSpecificities = new List<bool>();
            List<string> spliceJunctions = new List<string>();
            List<string> sortedBamFiles = new List<string>();
            List<string> dedupedBamFiles = new List<string>();

            // Trimming and strand specificity
            foreach (string[] fq in fastqs)
            {
                SkewerWrapper.Trim(bin, threads, 19, fq, out string[] trimmedFastqs, out string skewerLog);
                string[] fqForAlignment = trimmedFastqs;

                // Infer strand specificity
                bool localStrandSpecific = strandSpecific;
                if (inferStrandSpecificity)
                {
                    STARWrapper.SubsetFastqs(bin, fqForAlignment, readSubset, analysisDirectory, out string[] subsetFastqs);
                    if (useReadSubset) fqForAlignment = subsetFastqs;
                    string subsetOutPrefix = Path.Combine(Path.GetDirectoryName(subsetFastqs[0]), Path.GetFileNameWithoutExtension(subsetFastqs[0]));
                    WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignSubset.bash"), STARWrapper.BasicAlignReadCommands(bin, threads, genomeStarIndexDirectory, subsetFastqs, subsetOutPrefix, false, STARGenomeLoadOption.LoadAndKeep)).WaitForExit();
                    localStrandSpecific = RSeQCWrapper.CheckStrandSpecificity(bin, subsetOutPrefix + STARWrapper.BamFileSuffix, geneModelGtfOrGff, 0.8);
                }
                strandSpecificities.Add(localStrandSpecific);
                fastqsForAlignment.Add(fqForAlignment);
            }

            // Alignment
            List<string> alignmentCommands = new List<string>();
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SpliceJunctionFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.FirstPassAlignmentCommands(bin, threads, genomeStarIndexDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));
                }
                spliceJunctions.Add(outPrefix + STARWrapper.SpliceJunctionFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, genomeStarIndexDirectory));
            alignmentCommands.AddRange(STARWrapper.ProcessFirstPassSpliceCommands(spliceJunctions, out string spliceJunctionStartDatabase));
            int uniqueSuffix = 1;
            foreach (string f in fastqsForAlignment.SelectMany(f => f))
            {
                uniqueSuffix = uniqueSuffix ^ f.GetHashCode();
            }
            string secondPassGenomeDirectory = genomeStarIndexDirectory + "SecondPass" + uniqueSuffix.ToString();
            alignmentCommands.AddRange(STARWrapper.GenerateGenomeIndex(bin, threads, secondPassGenomeDirectory, new string[] { reorderedFasta }, geneModelGtfOrGff, spliceJunctionStartDatabase));
            foreach (string[] fq in fastqsForAlignment)
            {
                string outPrefix = Path.Combine(Path.GetDirectoryName(fq[0]), Path.GetFileNameWithoutExtension(fq[0]));
                if (!File.Exists(outPrefix + STARWrapper.SortedBamFileSuffix) || overwriteStarAlignment)
                {
                    alignmentCommands.AddRange(STARWrapper.AlignRNASeqReadsForVariantCalling(bin, threads, secondPassGenomeDirectory, fq, outPrefix, strandSpecificities[fastqsForAlignment.IndexOf(fq)], STARGenomeLoadOption.LoadAndKeep));

                }
                sortedBamFiles.Add(outPrefix + STARWrapper.SortedBamFileSuffix);
                dedupedBamFiles.Add(outPrefix + STARWrapper.DedupedBamFileSuffix);
            }
            alignmentCommands.AddRange(STARWrapper.RemoveGenome(bin, secondPassGenomeDirectory));
            WrapperUtility.GenerateAndRunScript(Path.Combine(bin, "scripts", "alignReads.bash"), alignmentCommands).WaitForExit();



        }

        public static void Test(string test)
        {
            string script_path = Path.Combine(test, "scripts", "test.sh");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(test),
                "echo " + "HaHa"
            }).WaitForExit();
        }
    }
}
