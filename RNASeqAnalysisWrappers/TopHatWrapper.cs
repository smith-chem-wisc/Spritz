using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNASeqAnalysisWrappers
{
    public class TopHatWrapper
    {
        public static void generate_index()
        {

        }

        public static void align()
        {

        }

        public static void install()
        {
            WrapperUtility.RunBashCommand("sudo", "apt-get install tophat");
        }
    }
}
