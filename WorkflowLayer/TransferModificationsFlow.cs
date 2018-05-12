using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    public class TransferModificationsFlow
    {
        public const string Command = "modtransfer";

        public List<string> TransferModifications(string spritzDirectory, string sourceXmlPath, List<string> destinationXmlPaths)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(spritzDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, new List<string>(), out Dictionary<string, Modification> un);
            List<string> outxmls = new List<string>();
            foreach (var xml in destinationXmlPaths)
            {
                string outxml = Path.Combine(Path.GetDirectoryName(xml), Path.GetFileNameWithoutExtension(xml) + ".withmods.xml");
                var newProts = ProteinAnnotation.TransferModifications(uniprot, ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, uniprotPtms, false, new List<string>(), out un));
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), newProts, outxml);
                outxmls.Add(outxml);
            }
            return outxmls;
        }

        public string TransferModifications(string spritzDirectory, string sourceXmlPath, string destinationXmlPath)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(spritzDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, new List<string>(), out Dictionary<string, Modification> un);
            string outxml = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".withmods.xml");
            var newProts = ProteinAnnotation.TransferModifications(uniprot, ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, new List<string>(), out un));
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), newProts, outxml);
            return outxml;
        }
    }
}