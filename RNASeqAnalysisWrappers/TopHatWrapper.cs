using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace RNASeqAnalysisWrappers
{
    public class TopHatWrapper
    {
        public void generate_index()
        {

        }

        public void align()
        {

        }

        public void install()
        {
            WrapperUtility.run_basic_command("sudo", "apt-get install tophat");
        }
    }
}
