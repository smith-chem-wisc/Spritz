﻿using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Web;
using UsefulProteomicsDatabases;

namespace SpritzModifications
{
    public static class ProteinAnnotation
    {
        /// <summary>
        /// Gets UniProt ptmlist
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static List<Modification> GetUniProtMods(string spritzDirectory)
        {
            Loaders.LoadElements();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(spritzDirectory, "PSI-MOD.obo.xml"));
            return Loaders.LoadUniprot(Path.Combine(spritzDirectory, "ptmlist.txt"), Loaders.GetFormalChargesDictionary(psiModDeserialized)).ToList();
        }

        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="proteogenomicProteins"></param>
        /// <param name="uniprotProteins"></param>
        /// <returns></returns>
        public static List<Protein> CombineAndAnnotateProteins(List<Protein> uniprotProteins, List<Protein> proteogenomicProteins)
        {
            var uniprotOrganism = uniprotProteins.Select(p => p.Organism).Distinct().ToList();
            List<Protein> newProteins = new();
            Dictionary<string, List<Protein>> dictUniprot = ProteinDictionary(uniprotProteins);
            Dictionary<string, List<Protein>> dictProteogenomic = ProteinDictionary(proteogenomicProteins);

            List<string> seqsInCommon = dictProteogenomic.Keys.Intersect(dictUniprot.Keys).ToList();
            List<string> pgOnlySeqs = dictProteogenomic.Keys.Except(dictUniprot.Keys).ToList();
            List<string> uniprotOnlySeqs = dictUniprot.Keys.Except(dictProteogenomic.Keys).ToList();

            foreach (string seq in seqsInCommon.Concat(pgOnlySeqs))
            {
                dictUniprot.TryGetValue(seq, out var uniprotProteinsSameSeq);
                dictProteogenomic.TryGetValue(seq, out var pgProteinsSameSeq);
                if (pgProteinsSameSeq == null) { continue; } // shouldn't happen

                var uniprot = CombineProteinsWithSameSeq(uniprotProteinsSameSeq); // may be null if a pg-only seq
                var pgProtein = CombineProteinsWithSameSeq(pgProteinsSameSeq);

                if (uniprot != null && uniprot.BaseSequence != pgProtein.BaseSequence) { throw new ArgumentException("Not all proteins have the same sequence"); }
                string proteinOrganism = null;
                if (uniprotOrganism.Count == 1)
                {
                    proteinOrganism = uniprotOrganism.FirstOrDefault();
                }
                else
                {
                    proteinOrganism = uniprot != null ? uniprot.Organism : pgProtein.Organism;
                }

                // keep all sequence variations
                newProteins.Add(new Protein(
                    pgProtein.BaseSequence,
                    uniprot != null ? uniprot.Accession : pgProtein.Accession, //comma-separated
                    organism: proteinOrganism,
                    name: pgProtein.Name, //comma-separated
                    fullName: uniprot != null ? uniprot.FullName : pgProtein.FullName, //comma-separated
                    isDecoy: pgProtein.IsDecoy,
                    isContaminant: pgProtein.IsContaminant,
                    sequenceVariations: pgProtein.SequenceVariations.OrderBy(v => v.OneBasedBeginPosition).ToList(),
                    spliceSites: pgProtein.SpliceSites.OrderBy(s => s.OneBasedBeginPosition).ToList(),

                    // combine these
                    geneNames: (uniprot != null ? uniprot.GeneNames : new List<Tuple<string, string>>()).Concat(pgProtein.GeneNames).ToList(),

                    // transfer these
                    oneBasedModifications: uniprot != null ? uniprot.OneBasedPossibleLocalizedModifications : new Dictionary<int, List<Modification>>(),
                    proteolysisProducts: uniprot != null ? uniprot.ProteolysisProducts.ToList() : new List<ProteolysisProduct>(),
                    databaseReferences: uniprot != null ? uniprot.DatabaseReferences.ToList() : new List<DatabaseReference>(),
                    disulfideBonds: uniprot != null ? uniprot.DisulfideBonds.ToList() : new List<DisulfideBond>()
                ));
            }

            // Some stats output
            var canonical = proteogenomicProteins.Select(p => p.NonVariantProtein).Distinct();
            Console.WriteLine($"{canonical.Count()}\tCanonincal proteins translated from gene model (without applied variations)");
            Console.WriteLine($"{seqsInCommon.Count}\tProteins with exact sequence match in UniProt");
            Console.WriteLine($"{pgOnlySeqs.Count}\tProteins without exact sequence match in UniProt");

            return newProteins;
        }

        /// <summary>
        /// Combine proteins containing the same sequence
        /// </summary>
        /// <param name="proteinsWithSameSequence"></param>
        /// <returns></returns>
        private static Protein CombineProteinsWithSameSeq(List<Protein> proteinsWithSameSequence)
        {
            if (proteinsWithSameSequence == null || proteinsWithSameSequence.Count <= 0)
            {
                return null; // no proteins to combine
            }

            var seq = new HashSet<string>(proteinsWithSameSequence.Select(p => p.BaseSequence));
            if (seq.Count > 1) { throw new ArgumentException("not all proteins have the same sequence"); }

            string accession = string.Join(",", proteinsWithSameSequence.Select(p => p.Accession));
            string fullName = string.Join(",", proteinsWithSameSequence.Select(p => p.FullName));
            var geneNames = new HashSet<Tuple<string, string>>(proteinsWithSameSequence.SelectMany(p => p.GeneNames)).ToList();
            if (accession.StartsWith("MSTRG")) // make it easier to interpret proteogenomic isoform IDs
            {
                string accession_pt1 = geneNames.FirstOrDefault(gn => gn.Item1 == "accession").Item2;
                string blastInformation = HttpUtility.UrlDecode(geneNames.FirstOrDefault(gn => gn.Item1 == "primary").Item2);
                blastInformation = HttpUtility.UrlDecode(geneNames.FirstOrDefault(gn => gn.Item1 == "primary").Item2);
                string[] blastResult = blastInformation.Split(",");
                string accession_pt2 = blastResult.Length > 2 ? $",{blastResult[2]}" : ",NA";
                accession = $"{accession_pt1}{accession_pt2}";
                fullName = $"{string.Join(",", proteinsWithSameSequence.Select(p => p.Accession))}|{blastInformation}";
            }

            return new Protein(
                seq.First(),
                accession,
                organism: proteinsWithSameSequence.First().Organism,
                name: string.Join(",", proteinsWithSameSequence.Select(p => p.Name)),
                fullName: fullName,
                isDecoy: proteinsWithSameSequence.All(p => p.IsDecoy),
                isContaminant: proteinsWithSameSequence.Any(p => p.IsContaminant),
                sequenceVariations: new HashSet<SequenceVariation>(proteinsWithSameSequence.SelectMany(p => p.SequenceVariations)).ToList(),
                spliceSites: new HashSet<SpliceSite>(proteinsWithSameSequence.SelectMany(p => p.SpliceSites)).ToList(),
                geneNames: geneNames,
                oneBasedModifications: CollapseMods(proteinsWithSameSequence),
                proteolysisProducts: new HashSet<ProteolysisProduct>(proteinsWithSameSequence.SelectMany(p => p.ProteolysisProducts)).ToList(),
                databaseReferences: new HashSet<DatabaseReference>(proteinsWithSameSequence.SelectMany(p => p.DatabaseReferences)).ToList(),
                disulfideBonds: new HashSet<DisulfideBond>(proteinsWithSameSequence.SelectMany(p => p.DisulfideBonds)).ToList()
            );
        }

        /// <summary>
        /// Collapses modification dictionaries from multiple proteins into one
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<int, List<Modification>> CollapseMods(IEnumerable<Protein> proteins)
        {
            Dictionary<int, HashSet<Modification>> result = new();
            if (proteins != null)
            {
                foreach (KeyValuePair<int, List<Modification>> kv in proteins.SelectMany(p => p.OneBasedPossibleLocalizedModifications).ToList())
                {
                    if (result.TryGetValue(kv.Key, out HashSet<Modification> modifications))
                    {
                        foreach (Modification mod in kv.Value)
                        {
                            modifications.Add(mod);
                        }
                    }
                    else
                    {
                        result[kv.Key] = new HashSet<Modification>(kv.Value);
                    }
                }
            }
            return result.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());
        }

        /// <summary>
        /// Creates a sequence-to-proteins dictionary
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<string, List<Protein>> ProteinDictionary(IEnumerable<Protein> proteins)
        {
            Dictionary<string, List<Protein>> dict = new();
            foreach (Protein p in proteins)
            {
                if (dict.TryGetValue(p.BaseSequence, out var prots)) { prots.Add(p); }
                else { dict[p.BaseSequence] = new List<Protein> { p }; }
            }
            return dict;
        }
    }
}