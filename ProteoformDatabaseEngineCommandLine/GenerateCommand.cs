using ManyConsole;

namespace ProteoformDatabaseEngine
{
    class GenerateCommand : ConsoleCommand
    {
        public GenerateCommand()
        {
            this.IsCommand("generate", "Generate a custom XML database.");

            // TODO: implement arguments using HasOption method from ManyConsole namespace
            // TODO: handle null input with ?? null-coalescing operator
            //HasOption("u|uniprot_xml=", "Enter the path of a uniprot database.", s => UniProtXml = s);
            //HasOption("e|ensembl_faa=", "Enter the path of an ensembl protein fasta database.", s => EnsemblFasta = s);
        }

        public override int Run(string[] remainingArguments)
        {

            return 0;
        }
    }
}
