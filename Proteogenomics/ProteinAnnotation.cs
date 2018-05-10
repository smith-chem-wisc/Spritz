using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class ProteinAnnotation
    {
        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="destination"></param>
        /// <param name="source"></param>
        /// <returns></returns>
        public static List<Protein> TransferModifications(List<Protein> source, List<Protein> destination)
        {
            List<Protein> newProteins = new List<Protein>();
            Dictionary<string, List<Protein>> dictSource = ProteinDictionary(source);
            Dictionary<string, List<Protein>> dictDestination = ProteinDictionary(destination);

            List<string> commonSeqs = dictDestination.Keys.Intersect(dictDestination.Keys).ToList();
            foreach (string seq in commonSeqs)
            {
                dictSource.TryGetValue(seq, out var sourceProteins);
                dictDestination.TryGetValue(seq, out var destinationProteins);
                if (sourceProteins == null || destinationProteins == null) { continue; }
                foreach (var destinationProtein in destinationProteins)
                {
                    newProteins.Add(
                        new Protein(

                            // keep these
                            seq,
                            destinationProtein.Accession,
                            name: destinationProtein.Name,
                            full_name: destinationProtein.FullName,
                            isDecoy: destinationProtein.IsDecoy,
                            isContaminant: destinationProtein.IsContaminant,
                            sequenceVariations: destinationProtein.SequenceVariations.ToList(),

                            // combine these
                            gene_names: destinationProtein.GeneNames.Concat(sourceProteins.SelectMany(p => p.GeneNames)).ToList(),

                            // transfer these
                            oneBasedModifications: CollapseMods(sourceProteins),
                            proteolysisProducts: new HashSet<ProteolysisProduct>(sourceProteins.SelectMany(p => p.ProteolysisProducts)).ToList(),
                            databaseReferences: new HashSet<DatabaseReference>(sourceProteins.SelectMany(p => p.DatabaseReferences)).ToList(),
                            disulfideBonds: new HashSet<DisulfideBond>(sourceProteins.SelectMany(p => p.DisulfideBonds)).ToList()
                        ));
                }
            }

            List<string> destinationOnly = dictDestination.Keys.Except(dictSource.Keys).ToList();
            List<string> sourceOnly = dictSource.Keys.Except(dictDestination.Keys).ToList();
            return newProteins;
        }

        /// <summary>
        /// Collapses modification dictionaries from multiple proteins into one
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<int, List<Modification>> CollapseMods(List<Protein> proteins)
        {
            Dictionary<int, HashSet<Modification>> result = new Dictionary<int, HashSet<Modification>>();
            foreach (var kv in proteins.SelectMany(p => p.OneBasedPossibleLocalizedModifications).ToList())
            {
                if (result.TryGetValue(kv.Key, out var modifications))
                {
                    foreach (var mod in kv.Value) { modifications.Add(mod); }
                }
                else
                {
                    result[kv.Key] = new HashSet<Modification>(kv.Value);
                }
            }
            return result.ToDictionary(kv => kv.Key, kv => kv.Value.ToList());
        }

        /// <summary>
        /// Creates a sequence-to-proteins dictionary
        /// </summary>
        /// <param name="proteins"></param>
        /// <returns></returns>
        private static Dictionary<string, List<Protein>> ProteinDictionary(List<Protein> proteins)
        {
            Dictionary<string, List<Protein>> dict = new Dictionary<string, List<Protein>>();
            foreach (Protein p in proteins)
            {
                if (dict.TryGetValue(p.BaseSequence, out var prots)) { prots.Add(p); }
                else { dict[p.BaseSequence] = new List<Protein> { p }; }
            }
            return dict;
        }
    }
}