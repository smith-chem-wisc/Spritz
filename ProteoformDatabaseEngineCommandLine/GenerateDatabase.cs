using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ManyConsole;

namespace ProteoformDatabaseEngine
{
    public class GenerateDatabaseCommand : ConsoleCommand
    {
        public GenerateDatabaseCommand()
        {
            this.IsCommand("generateDatabase", "Full run of alignment and RNA-Seq analysis to generate a proteomics database.");
        }

        public override int Run(string[] remainingArguments)
        {
            var commands = ProteoformDatabaseEngine.GetCommands();

            // USER PROMTS: commands
            List<string> commandList = new List<string> { "help", "generate", "gptmd", "exit" };
            string commandPrompt = String.Join(Environment.NewLine, new string[] { "Choose a command:",
                "help -- Display help menu for commandline arguments",
                "generate -- Generate a custom protein XML database",
                "gptmd -- Enter modifications predicted from first-pass GPTMD search results",
                "exit -- Exit this program" });

            // User prompts: parameters and input files
            List<string> trueFalseList = new List<string> { "T", "t", "True", "true", "TRUE", "F", "f", "False", "false", "FALSE" };
            List<string> trueList = new List<string> { "T", "t", "True", "true", "TRUE" };
            List<string> falseList = new List<string> { "F", "f", "False", "false", "FALSE" };
            string uniprotXmlPrompt = "Please enter the path for a UniProt XML.";
            //string ensemblPepAllPrompt = "Please enter the path for an Ensembl peptide fasta";
            string methionineOxidationPrompt = "";
            string carbamPrompt = "";
            string cleavedNTermMetPrompt = "";
            string lysineIsotopesPrompt = "";
            List<string> lysineOptList = new List<string>();
            string maxPtmsPrompt = "";
            int maxMaxPtms = 0;
            string proteoformXmlPrompt = "";
            string ensemblFaaPrompt = "";



            return 0;

        }
    }
}
