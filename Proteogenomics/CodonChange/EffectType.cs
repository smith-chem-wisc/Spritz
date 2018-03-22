namespace Proteogenomics
{
    public enum EffectType
    {
        // High impact
        // Order: Highest impact first
        CHROMOSOME_LARGE_DELETION //
        , CHROMOSOME_LARGE_INVERSION //
        , CHROMOSOME_LARGE_DUPLICATION //
        , GENE_REARRANGEMENT //
        , GENE_DELETED //
        , TRANSCRIPT_DELETED //
        , EXON_DELETED //
        , EXON_DELETED_PARTIAL //
        , GENE_FUSION //
        , GENE_FUSION_REVERESE //
        , GENE_FUSION_HALF //
        , FRAME_SHIFT //
        , STOP_GAINED //
        , STOP_LOST //
        , START_LOST //
        , SPLICE_SITE_ACCEPTOR //
        , SPLICE_SITE_DONOR //
        , RARE_AMINO_ACID //
        , EXON_DUPLICATION //
        , EXON_DUPLICATION_PARTIAL //
        , EXON_INVERSION //
        , EXON_INVERSION_PARTIAL //
        , PROTEIN_PROTEIN_INTERACTION_LOCUS //
        , PROTEIN_STRUCTURAL_INTERACTION_LOCUS //

        // Moderate impact
        // Order: Highest impact first
        // Note: Method Codon.effect() relies on this order for effect
        //       replacement (when 'allowReplace = true')
        , NON_SYNONYMOUS_CODING //
        , GENE_DUPLICATION //
        , TRANSCRIPT_DUPLICATION //
        , UTR_5_DELETED //
        , UTR_3_DELETED //
        , SPLICE_SITE_BRANCH_U12 //
        , GENE_INVERSION //
        , TRANSCRIPT_INVERSION //
        , CODON_INSERTION //
        , CODON_CHANGE_PLUS_CODON_INSERTION //
        , CODON_DELETION //
        , CODON_CHANGE_PLUS_CODON_DELETION //

        // Low impact
        // Order: Highest impact first
        , NON_SYNONYMOUS_STOP //
        , NON_SYNONYMOUS_START //
        , SPLICE_SITE_REGION //
        , SPLICE_SITE_BRANCH //
        , SYNONYMOUS_CODING //
        , SYNONYMOUS_START //
        , SYNONYMOUS_STOP //
        , CODON_CHANGE //
        , START_GAINED //
        , MOTIF //
        , MOTIF_DELETED //
        , FEATURE_FUSION //

        // Modifiers
        // Order: Highest impact first
        , UTR_5_PRIME //
        , UTR_3_PRIME //
        , REGULATION //
        , MICRO_RNA //
        , UPSTREAM //
        , DOWNSTREAM //
        , NEXT_PROT //
        , INTRON_CONSERVED //
        , INTRON //
        , INTRAGENIC //
        , INTERGENIC_CONSERVED //
        , INTERGENIC //
        , CDS //
        , EXON //
        , TRANSCRIPT //
        , GENE //
        , SEQUENCE //
        , CHROMOSOME_ELONGATION //
        , CUSTOM //
        , CHROMOSOME //
        , GENOME //
        , NONE //
    }
}