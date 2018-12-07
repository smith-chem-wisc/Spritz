using Proteomics;
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
                if (xml == null || !File.Exists(xml)) { continue; }
                string outxml = Path.Combine(Path.GetDirectoryName(xml), Path.GetFileNameWithoutExtension(xml) + ".withmods.xml");
                var nonVariantProts = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, uniprotPtms, false, null, out un).Select(p => p.NonVariantProtein).Distinct();
                var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, nonVariantProts.Concat(additionalProteins).ToList());
                ProteinDbWriter.WriteXmlDatabase(null, newProts, outxml);
                string outfasta = Path.Combine(Path.GetDirectoryName(xml), Path.GetFileNameWithoutExtension(xml) + ".spritz.fasta");
                ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).ToList(), outfasta, "|");
                outxmls.Add(outxml);
            }
            return outxmls;
        }

        public string TransferModifications(string spritzDirectory, string sourceXmlPath, string destinationXmlPath)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(spritzDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, null, out var un);
            string outxml = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".withmods.xml");
            var nonVariantProts = ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, null, out un).Select(p => p.NonVariantProtein).Distinct();
            var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, nonVariantProts.ToList());
            ProteinDbWriter.WriteXmlDatabase(null, newProts, outxml);
            string outfasta = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".spritz.fasta");
            var prot = newProts.FirstOrDefault(p => p.Accession.Contains("_"));
            ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).ToList(), outfasta, "|");
            return outxml;
        }
    }
}