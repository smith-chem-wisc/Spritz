using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class TranscriptQuantificationFlow
    {
        public string RsemReferenceIndexPrefix { get; private set; }
        public string RsemOutputPrefix { get; private set; }

        public void QuantifyTranscriptsFromSra(
            string binDirectory, string analysisDirectory, string referenceFastaPath, int threads, string geneModelPath, RSEMAlignerOption aligner, Strandedness strandedness,
            string sraAccession, bool doOutputBam)
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            RSEMWrapper rsem = new RSEMWrapper();
            sratoolkit.Fetch(binDirectory, sraAccession, analysisDirectory);
            string scriptName = Path.Combine(binDirectory, "scripts", "QuantifyTranscripts.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>(
                rsem.PrepareReferenceCommands(
                    binDirectory,
                    referenceFastaPath,
                    threads,
                    geneModelPath,
                    aligner)
                .Concat(
                rsem.CalculateExpressionCommands(
                    binDirectory,
                    RsemReferenceIndexPrefix,
                    threads,
                    aligner,
                    strandedness,
                    sratoolkit.FastqPaths,
                    doOutputBam))))
            .WaitForExit();

            RsemReferenceIndexPrefix = rsem.ReferencePrefix;
            RsemOutputPrefix = rsem.OutputPrefix;
        }

        public void QuantifyTranscripts(
            string binDirectory, string referenceFastaPath, int threads, string geneModelPath, RSEMAlignerOption aligner, Strandedness strandedness,
            string[] fastqPaths, bool doOutputBam)
        {
            SRAToolkitWrapper sratoolkit = new SRAToolkitWrapper();
            RSEMWrapper rsem = new RSEMWrapper();
            string scriptName = Path.Combine(binDirectory, "scripts", "QuantifyTranscripts.bash");
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>(
                rsem.PrepareReferenceCommands(
                    binDirectory,
                    referenceFastaPath,
                    threads,
                    geneModelPath,
                    aligner)
                .Concat(
                rsem.CalculateExpressionCommands(
                    binDirectory,
                    RsemReferenceIndexPrefix,
                    threads,
                    aligner,
                    strandedness,
                    fastqPaths,
                    doOutputBam))))
            .WaitForExit();

            RsemReferenceIndexPrefix = rsem.ReferencePrefix;
            RsemOutputPrefix = rsem.OutputPrefix;
        }
    }
}