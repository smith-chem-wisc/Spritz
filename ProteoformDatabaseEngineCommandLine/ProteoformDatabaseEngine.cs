using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using RNASeqAnalysisWrappers;
using CommandLine;

namespace ProteoformDatabaseEngineCommandLine
{
    class ProteoformDatabaseEngine
    {
        public static void Main(string[] args)
        {
            if (args.Contains("setup"))
            {
                WrapperUtility.Install(Environment.CurrentDirectory);
                return;
            }

            Options options = new Options();
            var isValid = Parser.Default.ParseArgumentsStrict(args, options);
            Fastq2ProteinsRunner.Run(options.Bin, options.GRCh37, options.GRCh38, options.Threads, new string[] { options.Fastq1, options.Fastq2 }, options.StrandSpecific, options.InferStrandSpecificity, options.GenomeStarIndexDirectory, options.GenomeFasta, options.GeneModelGtfOrGff, out string proteinDb);
            Console.WriteLine("ouptput database to " + proteinDb);
        }
    }
}
