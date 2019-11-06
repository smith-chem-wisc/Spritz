using Fclp;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace TransferUniProtModifications
{
    internal class TransferUniProtModifications
    {
        private static void Main(string[] args)
        {
            Console.WriteLine("Welcome to TransferModifications!");
            var p = new FluentCommandLineParser<ApplicationArguments>();

            p.Setup(arg => arg.UniProtXml)
                .As('x', "uniprot_xml")
                .Required()
                .WithDescription("UniProt protein XML file.");

            p.Setup(arg => arg.SpritzXml)
                .As('y', "spritz_xml")
                .Required()
                .WithDescription("Custom protein XML file, e.g. from Spritz.");

            p.SetupHelp("h", "help")
                .Callback(text => Console.WriteLine(text));

            var result = p.Parse(args);

            TransferModifications(p.Object.UniProtXml, p.Object.SpritzXml);
            DatabaseSummary(p.Object.UniProtXml, p.Object.SpritzXml);
        }

        public static string TransferModifications(string sourceXmlPath, string destinationXmlPath)
        {
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(Environment.CurrentDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, null, out var un);
            string outxml = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".withmods.xml");
            var nonVariantProts = ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, null, out un).Select(p => p.NonVariantProtein).Distinct();
            var newProts = ProteinAnnotation.CombineAndAnnotateProteins(uniprot, nonVariantProts.ToList());
            ProteinDbWriter.WriteXmlDatabase(null, newProts, outxml);
            string outfasta = Path.Combine(Path.GetDirectoryName(destinationXmlPath), Path.GetFileNameWithoutExtension(destinationXmlPath) + ".fasta");
            var prot = newProts.FirstOrDefault(p => p.Accession.Contains("_"));
            ProteinDbWriter.WriteFastaDatabase(newProts.SelectMany(p => p.GetVariantProteins()).ToList(), outfasta, "|");
            return outxml;
        }

        public static void DatabaseSummary(string sourceXmlPath, string destinationXmlPath)
        {
            var culture = CultureInfo.CurrentCulture;
            var uniprotPtms = ProteinAnnotation.GetUniProtMods(Environment.CurrentDirectory);
            var uniprot = ProteinDbLoader.LoadProteinXML(sourceXmlPath, true, DecoyType.None, uniprotPtms, false, null, out var un);
            var spritz = ProteinDbLoader.LoadProteinXML(destinationXmlPath, true, DecoyType.None, uniprotPtms, false, null, out un);
            var spritzCanonical = spritz.Select(p => p.NonVariantProtein).Distinct().ToList();
            int numberOfCanonicalProteinEntries = spritzCanonical.Count;
            int numberOfVariantProteinEntries = spritz.Count - spritzCanonical.Count;
            int synonymousCount = 0;
            int totalVariants = 0;
            int savCount = 0;
            int insertionCount = 0;
            int deletionCount = 0;
            int frameshiftCount = 0;
            int stopGainCount = 0;
            int stopLossCount = 0;
            int exonLossCount = 0;
            Dictionary<string, List<SequenceVariation>> allVariants = new Dictionary<string, List<SequenceVariation>>();
            foreach (var spritzEntry in spritz)
            {
                if (spritzEntry.AppliedSequenceVariations.Count() != 0)
                {
                    if (allVariants.ContainsKey(spritzEntry.NonVariantProtein.Accession))
                    {
                        foreach (var variant in spritzEntry.AppliedSequenceVariations)
                        {
                            if (!allVariants[spritzEntry.NonVariantProtein.Accession].Contains(variant))
                            {
                                allVariants[spritzEntry.NonVariantProtein.Accession].Add(variant);
                            }

                        }
                    }
                    else
                    {
                        allVariants.Add(spritzEntry.NonVariantProtein.Accession, spritzEntry.AppliedSequenceVariations);
                    }
                }
            }          
            foreach (var entry in allVariants)
            {
                foreach (var variant in entry.Value)
                {
                    if (culture.CompareInfo.IndexOf(variant.Description.Description, "synonymous_variant", CompareOptions.IgnoreCase) >= 0)
                    {
                        synonymousCount++;
                        totalVariants++;
                    }
                    if (culture.CompareInfo.IndexOf(variant.Description.Description, "missense_variant", CompareOptions.IgnoreCase) >= 0)
                    {
                        savCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "frameshift_variant", CompareOptions.IgnoreCase) >= 0)
                    {
                        frameshiftCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "stop_gained", CompareOptions.IgnoreCase) >= 0)
                    {
                        stopGainCount++;
                        totalVariants++;
                    }
                    else if ((culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_insertion", CompareOptions.IgnoreCase) >= 0) || (culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_insertion", CompareOptions.IgnoreCase) >= 0))
                    {
                        insertionCount++;
                        totalVariants++;
                    }
                    else if ((culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_deletion", CompareOptions.IgnoreCase) >= 0) || (culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_deletion", CompareOptions.IgnoreCase) >= 0))
                    {
                        deletionCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "exon_loss_variant", CompareOptions.IgnoreCase) >= 0)
                    {
                        exonLossCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "stop_loss", CompareOptions.IgnoreCase) >= 0)
                    {
                        stopLossCount++;
                        totalVariants++;
                    }
                }

            }

            string[] summary = new string[20];
            summary[0] = "Spritz Database Summary";
            summary[1] = "--------------------------------------------------------------";
            summary[2] = "Total number of protein entries in the database: " + spritz.Count();
            summary[3] = "Total number of canonical protein entries in the database: " + numberOfCanonicalProteinEntries;
            summary[4] = "Total number of variant containing protein entries in the database: " + numberOfVariantProteinEntries;
            summary[5] = "  Total number of unique variants in the database: " + totalVariants;
            summary[6] = "      Total number of  unique synonymous variants in the database: " + synonymousCount;
            summary[7] = "      Total number of unique nonsynonymous variants in the database: " + (totalVariants - synonymousCount);
            summary[8] = "          Number of unique SAVs in the database: " + savCount;
            summary[9] = "          Number of unique frameshift variants in the database: " + frameshiftCount;
            summary[10] = "         Number of unique insertion variants in the database: " + insertionCount;
            summary[11] = "         Number of unique deletion variants in the database: " + deletionCount;
            summary[12] = "         Number of unique stop gain variants in the database: " + stopGainCount;
            summary[13] = "         Number of unique stop loss variants in the database: " + stopLossCount;
            summary[14] = "         Number of unique exon loss variants int he database: " + exonLossCount;

            File.WriteAllLines(Path.Combine(Path.GetDirectoryName(destinationXmlPath), "SpritzDatabaseSummary.txt"), summary);

        }

        public class ApplicationArguments
        {
            public string SpritzXml { get; set; }
            public string ReferenceGeneModel { get; set; }
            public string UniProtXml { get; set; }
        }
    }
}