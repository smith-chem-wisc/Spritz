using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using CommandLine;

namespace ProteoformDatabaseEngineCommandLine
{
    class Options
    {
        [Option('c', "command", Required = true, HelpText = "command: (1) setup or (2) run")]
        public string Command { get; set; }

        [Option('b', "binDirectory", Required =true, HelpText ="bin directory for PDE")]
        public string Bin { get; set; }

        [Option('7', "GRCh37", Required = false, HelpText = "use GRCH37 references", DefaultValue = false)]
        public bool GRCh37 { get; set; }

        [Option('8', "GRCh38", Required = false, HelpText = "use GRCH38 references", DefaultValue = true)]
        public bool GRCh38 { get; set; }

        [Option('t', "threads", Required = false, HelpText = "number of threads to use")]
        public int Threads { get; set; } = Environment.ProcessorCount;

        [Option('1', "fq1", Required = true, HelpText = "fastq for single-end or for pair1")]
        public string Fastq1 { get; set; }

        [Option('2', "fq2", Required = false, HelpText = "fastq pair2")]
        public string Fastq2 { get; set; }

        [Option('s', "strandSpecific", Required = false, HelpText = "strandedness of the protocol", DefaultValue =false)]
        public bool StrandSpecific { get; set; }

        [Option('i', "inferStrandedness", Required = true, HelpText = "infer the strandedness with a sample of the fastq files", DefaultValue =true)]
        public bool InferStrandSpecificity { get; set; }

        [Option("genomeDir", Required = true, HelpText = "STAR genome directory", DefaultValue = true)]
        public string GenomeStarIndexDirectory { get; set; }

        [Option("genomeFasta", Required = true, HelpText = "genomeFasta", DefaultValue = true)]
        public string GenomeFasta { get; set; }

        [Option("geneModelGtfOrGff", Required = true, HelpText = "geneModelGtfOrGff", DefaultValue = true)]
        public string GeneModelGtfOrGff { get; set; }
    }
}
