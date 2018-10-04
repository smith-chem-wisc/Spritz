using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    public class TransferModificationsFlow
    {
        public const string Command = "modtransfer";

        public List<string> TransferModifications(string spritzDirectory, string sourceXmlPath, List<string> destinationXmlPaths, List<Protein> additionalProteins)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(spritzDirectory);
            List<string> outxmls = new List<string>();

            var uniprot = File.Exists(sourceXmlPath) ?
                ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, null, out Dictionary<string, Modification> un) :
                new List<Protein>();

            foreach (var xml in destinationXmlPaths)
            {
                string outxml = Path.Combine(Path.GetDirectoryName(xml), Path.GetFileNameWithoutExtension(xml) + ".withmods.xml");
                var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, uniprotPtms, false, null, out un).Concat(additionalProteins).ToList());
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), newProts, outxml);
                string outfasta = Path.Combine(Path.GetDirectoryName(xml), Path.GetFileNameWithoutExtension(xml) + ".spritz.fasta");
                ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).OfType<Protein>().ToList(), outfasta, "|");
                outxmls.Add(outxml);
            }
            return outxmls;
        }

        public string TransferModifications(string spritzDirectory, string sourceXmlPath, string destinationXmlPath)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(spritzDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, new List<string>(), out Dictionary<string, Modification> un);
            string outxml = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".withmods.xml");
            var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, new List<string>(), out un));
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), newProts, outxml);
            string outfasta = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".spritz.fasta");
            ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).OfType<Protein>().ToList(), outfasta, "|");
            return outxml;
        }
    }
}