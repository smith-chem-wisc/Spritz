using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ToolWrapperLayer
{
    /// <summary>
    /// RSEM is a program for calculating RNA-Seq Abundancy by Estimation Maximization.
    /// </summary>
    public class RSEMWrapper :
        IInstallable
    {

        #region Installation Methods

        /// <summary>
        /// Writes an installation script for RSEM.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string binDirectory)
        {
            return null;
        }

        /// <summary>
        /// Writes a script for removing RSEM.
        /// </summary>
        /// <param name="binDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string binDirectory)
        {
            return null;
        }

        #endregion Installation Methods

    }
}
