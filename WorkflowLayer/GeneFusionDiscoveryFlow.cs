using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class GeneFusionDiscoveryFlow
     : SpritzFlow
    {
        public const string Command = "fusion";

        public GeneFusionDiscoveryFlow()
            : base(MyWorkflow.TranscriptQuantification)
        {
        }

        public GeneFusionDiscoveryParameters Parameters { get; set; } = new GeneFusionDiscoveryParameters();
        public List<Protein> FusionProteins { get; set; } = new List<Protein>();

        public void DiscoverGeneFusions()
        {
            HashSet<string> usedFusionProteinAccessions = new HashSet<string>();
            foreach (string[] fastqs in Parameters.Fastqs)
            {
                // Run workflow
                STARFusionWrapper fusion = new STARFusionWrapper();
                string scriptName = WrapperUtility.GetAnalysisScriptPath(Parameters.AnalysisDirectory, "GeneFusionWorkflow.bash");
                var referenceCommands = fusion.DownloadPrecompiledReference(Parameters.SpritzDirectory, Parameters.Reference);
                var calculateCommands = fusion.RunStarFusion(Parameters.SpritzDirectory, Parameters.AnalysisDirectory, Parameters.Threads, fastqs);
                WrapperUtility.GenerateAndRunScript(scriptName, new List<string>(referenceCommands.Concat(calculateCommands))).WaitForExit();

                // Process results
                FusionProteins.AddRange(fusion.ParseCodingEffect(Path.Combine(fusion.OutputDirectoryPath, fusion.CodingEffectFilename),
                    Parameters.MinPeptideLength, Parameters.Organism, usedFusionProteinAccessions));
            }
        }

        /// <summary>
        /// Run this workflow (for GUI)
        /// </summary>
        /// <param name="parameters"></param>
        protected override void RunSpecific(ISpritzParameters parameters)
        {
            Parameters = (GeneFusionDiscoveryParameters)parameters;
            DiscoverGeneFusions();
        }
    }
}