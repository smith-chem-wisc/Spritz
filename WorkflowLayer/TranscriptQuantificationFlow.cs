using System.Collections.Generic;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class TranscriptQuantificationFlow
        : SpritzFlow
    {
        public const string Command = "quantify";

        public TranscriptQuantificationFlow()
            : base(MyWorkflow.TranscriptQuantification)
        {
        }

        public TranscriptQuantificationParameters Parameters { get; set; }
        public string RsemReferenceIndexPrefix { get; private set; }
        public string RsemOutputPrefix { get; private set; }

        public void QuantifyTranscripts()
        {
            RSEMWrapper rsem = new RSEMWrapper();
            string scriptName = WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "QuantifyTranscripts.bash");
            var referenceCommands = rsem.PrepareReferenceCommands(
                    Parameters.SpritzDirectory,
                    Parameters.ReferenceFastaPath,
                    Parameters.Threads,
                    Parameters.GeneModelPath,
                    Parameters.Aligner);
            var calculateCommands = rsem.CalculateExpressionCommands(
                    Parameters.SpritzDirectory,
                    rsem.ReferenceIndexPrefix,
                    Parameters.Threads,
                    Parameters.Aligner,
                    Parameters.Strandedness,
                    Parameters.Fastq,
                    Parameters.DoOutputQuantificationBam);
            WrapperUtility.GenerateAndRunScript(scriptName, new List<string>(referenceCommands.Concat(calculateCommands))).WaitForExit();

            RsemReferenceIndexPrefix = rsem.ReferenceIndexPrefix;
            RsemOutputPrefix = rsem.OutputPrefix;

            RsemReferenceIndexPrefix = rsem.ReferenceIndexPrefix;
            RsemOutputPrefix = rsem.OutputPrefix;
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (TranscriptQuantificationParameters)parameters;
            QuantifyTranscripts();
        }
    }
}