﻿namespace SpritzBackend
{
    public static class SpritzCmdAppArgInfoStrings
    {
        public static readonly char AnalysisDirectoryShort = 'a';
        public static readonly string AnalysisDirectoryLong = "analysisDirectory";
        public static readonly string AnalysisDirectoryDesc =
            "Directory in which the results will be stored.";

        public static readonly char AnalyzeVariantsShort = 'b';
        public static readonly string AnalyzeVariantsLong = "analyzeVariants";
        public static readonly string AnalyzeVariantsDesc =
            "Analyze protein coding variations and include them in the proteogenomic database. " +
            "Note: Specifying false for both AnalyzeVariants and AnalyzeIsoforms will generate a " +
            "reference proteogenomic database from the Ensembl references without any variants.";

        public static readonly char AnalyzeIsoformsShort = 'c';
        public static readonly string AnalyzeIsoformsLong = "analyzeIsoforms";
        public static readonly string AnalyzeIsoformsDesc =
            "Analyze alternative splicing events and include them in the proteogenomic database.";

        public static readonly char QuantifyShort = 'd';
        public static readonly string QuantifyLong = "doQuantification";
        public static readonly string QuantifyDesc =
            "Quantify genes and isoforms after constructing proteogenomic database.";

        public static readonly char AvailableReferencesShort = 'x';
        public static readonly string AvailableReferencesLong = "availableReferences";
        public static readonly string AvailableReferencesDesc =
            "Save a comma-separated file with available references to analysis directory. Then, exit.";

        public static readonly char AnalysisSetupShort = 'y';
        public static readonly string AnalysisSetupLong = "analysisSetup";
        public static readonly string AnalysisSetupDesc =
            "Perform setup required for protected access servers, given reference, species, and organism specified.";

        public static readonly char Fastq1Short = 'i';
        public static readonly string Fastq1Long = "fastq1";
        public static readonly string Fastq1Desc = "Comma-separated list of paths to first mate pair fastq files.";

        public static readonly char Fastq2Short = 'j';
        public static readonly string Fastq2Long = "fastq2";
        public static readonly string Fastq2Desc =
            "Comma-separated list of paths to second mate pair fastq files.";

        public static readonly char Fastq1SingleEndShort = 'f';
        public static readonly string Fastq1SingleEndLong = "fastq1SingleEnd";
        public static readonly string Fastq1SingleEndDesc =
            "Comma-separated list of paths to single-end fastq files.";

        public static readonly char SraAccessionShort = 's';
        public static readonly string SraAccessionLong = "sraAccessionPairedEnd";
        public static readonly string SraAccessionDesc =
            "Comma-separated list of SRA accessions for paired-end experiments to download.";

        public static readonly char SraAccessionSingleEndShort = 't';
        public static readonly string SraAccessionSingleEndLong = "sraAccessionSingleEnd";
        public static readonly string SraAccessionSingleEndDesc =
            "Comma-separated list of SRA accessions for single-end experiments to download.";

        public static readonly char ThreadsShort = 'p';
        public static readonly string ThreadsLong = "threads";
        public static readonly string ThreadsDesc =
            "Number of processors to use for analysis.";

        public static readonly char ReferenceShort = 'r';
        public static readonly string ReferenceLong = "reference";
        public static readonly string ReferenceDesc =
            "Reference to use, e.g. release-96,homo_sapiens,human,GRCh38. " +
            "Copy-paste a line from the file you get with the -x option that retrieves available references.";
    }
}