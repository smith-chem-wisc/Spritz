using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ToolWrapperLayer
{
    /// <summary>
    /// MELTING USAGE
    /// 
    ///    Mandatory Parameters:
    ///    
    ///    -S [XXXXXXXXXX] Nucleic acid sequence
    ///    
    ///    -C[XXXXXXXXXX] Complementary sequence, mandatory if mismaches
    ///    
    ///    -H[xxxxxx] Type of hybridisation (exemple dnadna)
    ///    
    ///    -N[x.xe - x] Sodium concentration in mol.l-1.
    ///    
    ///    -k[x.xe - x] Potassium concentration in mol.l-1.
    ///    
    ///    -t[x.xe - x] Tris concentration in mol.l-1. The Tri+ concentration is about half of total Tris concentration
    ///        
    ///    -G[x.xe - x] Magnesium concentration in mol.l-1.
    ///
    ///    -----------------------------------------------------------
    ///    
    ///    Optional Parameters:
    ///    
    ///    -A[xxxxxx.nn] Name of a file containing alternative nn parameters
    ///                   Defaults are: DNA/DNA: all97a.nn
    ///                                  DNA/RNA: sug95a.nn
    ///                                  RNA/RNA: xia98a.nn
    ///                                  
    ///    -D[xxxxxx.nn] Name of a file containing nn parameters for dangling ends
    ///                    Default is dnadnade.nn
    ///                    
    ///     
    ///    -F[x.xx] Correction for the concentration of nucleic acid
    ///                    Default is DEFAULT_NUC_CORR
    ///                    
    ///    -h Displays this help and quit
    ///     
    ///    
    ///    -I[XXXXXX] Name of an input file setting up the options
    ///     
    ///    -K Salt correction.Default is san98a
    ///     
    ///    -L Displays legal information and quit
    ///    
    ///    -M[xxxxxx.nn] Name of a file containing nn parameters for mismatches
    ///        Default is dnadnamm.nn
    ///        
    ///    -i[xxxxxx.nn] Name of a file containing nn parameters for inosine mismatches
    ///        Defaults are: DNA/DNA: san05a.nn
    ///                      DNA/RNA: san05a.nn
    ///                      RNA/RNA: bre07a.nn
    ///    
    ///    -O[XXXXXX] Name of an output file (the name can be omitted)
    ///    
    ///    -P [x.xe-x] Concentration of single strand nucleic acid in mol.l-1. Mandatory
    ///    
    ///    -p Return path where to find the calorimetric tables
    ///     
    ///    -q Quiet. Switch off interactive correction of parameters
    ///    
    ///    -T [XXX] Threshold for approximative computation
    ///    
    ///    -v Switch ON the verbose mode, issuing lot more info
    ///                    (if already ON, switch if OFF). Default is OFF
    ///    
    ///    -V Print the version number
    ///    
    ///    -x Force to compute an approximative tm
    ///    
    ///    More information is available in the user-guide. Type `man melting'
    ///     to access it, or consult one of the melting.xxx files, where xxx
    ///     states for lat1 (isolatin1 text), ps (postscript), pdf or html.
    /// </summary>
    public class MeltingWrapper
    {

    }
}
