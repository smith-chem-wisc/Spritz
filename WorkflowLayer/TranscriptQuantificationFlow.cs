using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class TranscriptQuantificationFlow
    {
        public static void QuantifyTranscripts(
            string binDirectory, string referenceFastaPath, int threads, string geneModelPath, RSEMAlignerOption aligner, Strandedness strandedness,
            string[] fastqPaths, bool doOutputBam, out string rsemReferencePrefix, out string outputPrefix)
        {
            WrapperUtility.GenerateAndRunScript(Path.Combine(binDirectory, "scripts", "QuantifyTranscripts.bash"), new List<string>(
                RSEMWrapper.PrepareReferenceCommands(
                    binDirectory,
                    referenceFastaPath,
                    threads,
                    geneModelPath,
                    aligner,
                    out rsemReferencePrefix)
                .Concat(
                RSEMWrapper.CalculateExpressionCommands(
                    binDirectory,
                    rsemReferencePrefix,
                    threads,
                    aligner,
                    strandedness,
                    fastqPaths,
                    doOutputBam,
                    out outputPrefix))))
            .WaitForExit();
        }

        public static void QuantifyTranscriptsFromSra(
            string binDirectory, string analysisDirectory, string referenceFastaPath, int threads, string geneModelPath, RSEMAlignerOption aligner, Strandedness strandedness,
            string sraAccession, bool doOutputBam, out string rsemReferencePrefix, out string outputPrefix)
        {
            SRAToolkitWrapper.Fetch(binDirectory, sraAccession, analysisDirectory, out string[] fastqPaths, out string logPath);
            WrapperUtility.GenerateAndRunScript(Path.Combine(binDirectory, "scripts", "QuantifyTranscripts.bash"), new List<string>(
                RSEMWrapper.PrepareReferenceCommands(
                    binDirectory,
                    referenceFastaPath,
                    threads,
                    geneModelPath,
                    aligner,
                    out rsemReferencePrefix)
                .Concat(
                RSEMWrapper.CalculateExpressionCommands(
                    binDirectory,
                    rsemReferencePrefix,
                    threads,
                    aligner,
                    strandedness,
                    fastqPaths,
                    doOutputBam,
                    out outputPrefix))))
            .WaitForExit();
        }
    }
}