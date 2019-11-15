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
                .WithDescription("Custom protein XML file, e.g. from Spritz.");

            p.Setup(arg => arg.FusionCodingEffects)
                .As('f', "fusion_coding_effect")
                .WithDescription("Coding effects from STAR-Fusion, comma separated");

            p.SetupHelp("h", "help")
                .Callback(text => Console.WriteLine(text));

            var result = p.Parse(args);

            TransferModifications(p.Object.UniProtXml, p.Object.SpritzXml ?? ProteinAnnotation.ParseCodingEffectsToXml(p.Object.FusionCodingEffects));
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
            int missenseSnvCount = 0;
            int missenseMnvCount = 0;
            int insertionCount = 0;
            int deletionCount = 0;
            int frameshiftCount = 0;
            int stopGainCount = 0;
            int stopLossCount = 0;
            Dictionary<string, List<SequenceVariation>> allVariants = new Dictionary<string, List<SequenceVariation>>();
            foreach (var spritzEntry in spritz)
            {
                if (spritzEntry.AppliedSequenceVariations.Count != 0)
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
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "missense_variant", CompareOptions.IgnoreCase) >= 0 &&
                        variant.Description.ReferenceAlleleString.Length == 1 && variant.Description.AlternateAlleleString.Length == 1)
                    {
                        missenseSnvCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "missense_variant", CompareOptions.IgnoreCase) >= 0)
                    {
                        missenseMnvCount++;
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
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_insertion", CompareOptions.IgnoreCase) >= 0 || culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_insertion", CompareOptions.IgnoreCase) >= 0)
                    {
                        insertionCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "conservative_inframe_deletion", CompareOptions.IgnoreCase) >= 0 || culture.CompareInfo.IndexOf(variant.Description.Description, "disruptive_inframe_deletion", CompareOptions.IgnoreCase) >= 0)
                    {
                        deletionCount++;
                        totalVariants++;
                    }
                    else if (culture.CompareInfo.IndexOf(variant.Description.Description, "stop_loss", CompareOptions.IgnoreCase) >= 0)
                    {
                        stopLossCount++;
                        totalVariants++;
                    }
                }
            }

            Console.WriteLine($"Spritz Database Summary");
            Console.WriteLine($"--------------------------------------------------------------");
            Console.WriteLine($"{numberOfCanonicalProteinEntries}\tTotal number of canonical protein entries (before applying variations)");
            Console.WriteLine($"{spritz.Count}\tTotal number of protein entries");
            Console.WriteLine($"{numberOfVariantProteinEntries}\tTotal number of variant containing protein entries");
            Console.WriteLine($"{totalVariants}\tTotal number of unique variants");
            Console.WriteLine($"{synonymousCount}\tTotal number of unique synonymous variants");
            Console.WriteLine($"{(totalVariants - synonymousCount)}\tTotal number of unique nonsynonymous variants");
            Console.WriteLine($"{missenseSnvCount}\tNumber of unique SNV missense variants");
            Console.WriteLine($"{missenseMnvCount}\tNumber of unique MNV missense variants");
            Console.WriteLine($"{frameshiftCount}\tNumber of unique frameshift variants");
            Console.WriteLine($"{insertionCount}\tNumber of unique insertion variants");
            Console.WriteLine($"{deletionCount}\tNumber of unique deletion variants");
            Console.WriteLine($"{stopGainCount}\tNumber of unique stop gain variants");
            Console.WriteLine($"{stopLossCount}\tNumber of unique stop loss variants");
        }

        public class ApplicationArguments
        {
            public string SpritzXml { get; set; }
            public string ReferenceGeneModel { get; set; }
            public string UniProtXml { get; set; }
            public string FusionCodingEffects { get; set; }
        }
    }
}