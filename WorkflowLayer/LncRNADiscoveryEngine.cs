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
    public class LncRNADiscoveryEngine
    {
        
        public static void Test(string test)
        {
            string script_path = Path.Combine(test, "scripts", "test.sh");
            WrapperUtility.GenerateAndRunScript(script_path, new List<string>
            {
                "cd " + WrapperUtility.ConvertWindowsPath(test),
                "echo " + "HaHa"
            }).WaitForExit();
        }

        public static void RunLncRNADiscoveryFromFastq(string bin, string analysisDirectory, int threads, List<string[]> fastqs, bool strandSpecific, bool inferStrandSpecificity, bool overwriteStarAlignment, string genomeStarIndexDirectory, string genomeFasta, string geneModelGtfOrGff, bool useReadSubset = false, int readSubset = 300000)
        {

        }

        public static void RunCufflink()
        {

        }

        public static void Rungtf2bed()
        {

        }

        public static void RunEnsembl2UCSC()
        {

        }

        public static void RunSlncky()
        {

        }

    }
}
