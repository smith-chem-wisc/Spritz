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
    public class LncRNADiscoveryFlow : SpritzFlow
    {
        public LncRNADiscoveryFlow() : base(MyWorkflow.LncRnaDiscovery)
        {
            Parameters = new Parameters();
        }   

        public Parameters Parameters { get; set; }
       
        public static void Test(string test)
        {
            string script_path = Path.Combine(test, "scripts", "test.sh");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(test),
                "echo " + "HaHa"
            }).WaitForExit();
        }

        protected override void RunSpecific(string OutputFolder, List<string> genomeFastaList, List<string> geneSetList, List<string> rnaSeqFastqList)
        {
            Test(OutputFolder);
        }
    }
}
