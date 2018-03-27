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
            Dictionary<string, Protein> dictDestination = destination.ToDictionary(p => p.BaseSequence, p => p);
            Dictionary<string, Protein> dictSource = source.ToDictionary(p => p.BaseSequence, p => p);

            List<string> commonSeqs = dictDestination.Keys.Intersect(dictDestination.Keys).ToList();
            foreach (string seq in commonSeqs)
            {
                Protein source_protein = dictSource[seq];
                Protein destination_protein = dictDestination[seq];
                newProteins.Add(
                    new Protein(
                        seq,
                        destination_protein.Accession,
                        gene_names: source_protein.GeneNames.ToList(),
                        oneBasedModifications: source_protein.OneBasedPossibleLocalizedModifications,
                        proteolysisProducts: source_protein.ProteolysisProducts.ToList(),
                        name: destination_protein.Name,
                        full_name: destination_protein.FullName,
                        isDecoy: destination_protein.IsDecoy,
                        isContaminant: destination_protein.IsContaminant,
                        databaseReferences: source_protein.DatabaseReferences.ToList(),
                        sequenceVariations: source_protein.SequenceVariations.ToList()));
            }

            List<string> destinationOnly = dictDestination.Keys.Except(dictSource.Keys).ToList();
            List<string> sourceOnly = dictSource.Keys.Except(dictDestination.Keys).ToList();
            return newProteins;
        }
    }
}