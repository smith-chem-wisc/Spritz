using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ManyConsole;

namespace ProteoformDatabaseEngine
{
    class ProteoformDatabaseEngine
    {
        public static void Main(string[] args)
        {
            var commands = GetCommands();
            if (args.Length > 0)
            {
                ConsoleCommandDispatcher.DispatchCommand(commands, args, Console.Out);
            }
            else
            {
                ConsoleCommandDispatcher.DispatchCommand(commands, new string[] { "user_interface" }, Console.Out);
            }
        }

        public static IEnumerable<ConsoleCommand> GetCommands()
        {
            return ConsoleCommandDispatcher.FindCommandsInSameAssemblyAs(typeof(ProteoformDatabaseEngine));
        }
    }
}
