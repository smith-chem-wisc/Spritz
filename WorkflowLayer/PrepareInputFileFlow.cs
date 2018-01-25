using System;
using System.IO;
using Proteogenomics;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class PrepareInputFileFlow
    {
        public static void PrepareEnsemblGenomeFasta(string genomeFasta, out Genome ensemblGenome, out string reorderedFasta)
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
            reorderedFasta = Path.Combine(Path.GetDirectoryName(genomeFasta), Path.GetFileNameWithoutExtension(genomeFasta) + ".karyotypic.fa");
            ensemblGenome = new Genome(genomeFasta);
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
        }
    }
}
