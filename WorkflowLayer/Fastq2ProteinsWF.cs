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
    public class Fastq2ProteinsWF : SpritzWorkflow
    {
        public Fastq2ProteinsWF() : base(MyWorkflow.Fastaq2Proteins)
        {
            rnaSeqAlignParameters = new RnaSeqAlignParameters();
        }

        public RnaSeqAlignParameters rnaSeqAlignParameters { get; set; }

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