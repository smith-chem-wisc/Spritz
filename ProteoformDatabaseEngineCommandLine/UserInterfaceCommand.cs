using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ManyConsole;

namespace ProteoformDatabaseEngine
{
    class UserInterfaceCommand : ConsoleCommand
    {
        public UserInterfaceCommand()
        {
            this.IsCommand("user_interface", "Guided console command setup.");
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

            // RUN COMMANDS after prompting and parsing arguments
            string command = GetUserResponse(commandPrompt);
            if (command == "exit") { }
            else if (command == "help") { ConsoleCommandDispatcher.DispatchCommand(commands, remainingArguments, Console.Out); }
            else
            {
                string[] newargs = new string[0];

                if (command == "generate")
                {
                    bool methionineOxidation = trueList.Contains(GetUserResponse(methionineOxidationPrompt, trueFalseList));
                    bool carbam = trueList.Contains(GetUserResponse(carbamPrompt, trueFalseList));
                    bool cleaved_met = trueList.Contains(GetUserResponse(cleavedNTermMetPrompt, trueFalseList));
                    string lysine_isotopes = GetUserResponse(lysineIsotopesPrompt, lysineOptList);
                    int maxPtms = GetUserResponse(maxPtmsPrompt, Enumerable.Range(0, maxMaxPtms).ToList());

                    string proteoformXml = GetUserResponse(proteoformXmlPrompt);
                    string uniprotXml = GetUserResponse(uniprotXmlPrompt);
                    string ensemblFaa = GetUserResponse(ensemblFaaPrompt);
                    newargs = new string[7] { command, "-f", proteoformXml, "-u", uniprotXml, "-e", ensemblFaa }; //TODO: test that this works to specify the correct args for ManyConsole
                }

                ConsoleCommandDispatcher.DispatchCommand(commands, newargs, Console.Out); // TODO: implement error notification if empty newargs at this point
            }

            return 0;

        }

        //User response methods (+2 overload)
        public string GetUserResponse(string prompt, List<string> acceptableAnswers)
        {
            string response = null;
            do
            {
                Console.WriteLine(prompt);
                response = Console.ReadLine().Trim();
            } while (!acceptableAnswers.Contains(response));
            return response;
        }
        public int GetUserResponse(string prompt, List<int> acceptableAnswers)
        {
            int response = Int32.MinValue;
            do
            {
                try
                {
                    Console.WriteLine(prompt);
                    response = Int32.Parse(Console.ReadLine().Trim());
                }
                catch { }
            } while (!acceptableAnswers.Contains(response));
            return response;
        }
        public string GetUserResponse(string prompt)
        {
            Console.WriteLine(prompt);
            return Console.ReadLine().Trim();
        }
    }
}
