namespace Proteogenomics
{
    public enum ErrorWarningType
    {
        INFO_REALIGN_3_PRIME // 		Variant has been realigned to the most 3-prime position within the transcript. This is usually done to to comply with HGVS specification to always report the most 3-prime annotation.
        , WARNING_SEQUENCE_NOT_AVAILABLE // The exon does not have reference sequence information
        , WARNING_REF_DOES_NOT_MATCH_GENOME // Sequence reference does not match variant's reference (alignment problem?)
        , WARNING_TRANSCRIPT_INCOMPLETE // Number of coding bases is NOT multiple of 3, so there is missing information for at least one codon.
        , WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS // Multiple STOP codons found in a CDS (usually indicates frame errors un one or more exons)
        , WARNING_TRANSCRIPT_NO_START_CODON // Start codon does not match any 'start' codon in the CodonTable
        , WARNING_TRANSCRIPT_NO_STOP_CODON // Stop codon does not match any 'stop' codon in the CodonTable
        , ERROR_CHROMOSOME_NOT_FOUND // Chromosome name not found. Typically due to mismatch in chromosome naming conventions between variants file and database, but can be a more severa problem (different reference genome)
        , ERROR_OUT_OF_CHROMOSOME_RANGE // Variant is outside chromosome
        , ERROR_OUT_OF_EXON //
        , ERROR_MISSING_CDS_SEQUENCE // Missing coding sequence information
    }
}