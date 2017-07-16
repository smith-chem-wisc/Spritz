using Proteomics;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class SequenceSimilarity
    {
        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="destination"></param>
        /// <param name="source"></param>
        /// <returns></returns>
        public static List<Protein> TransferModifications(List<Protein> source, List<Protein> destination)
        {
            List<Protein> new_proteins = new List<Protein>();
            Dictionary<string, Protein> dict_destination = destination.ToDictionary(p => p.BaseSequence, p => p);
            Dictionary<string, Protein> dict_source = source.ToDictionary(p => p.BaseSequence, p => p);

            List<string> common_seqs = dict_destination.Keys.Intersect(dict_destination.Keys).ToList();
            foreach (string seq in common_seqs)
            {
                Protein source_protein = dict_source[seq];
                Protein destination_protein = dict_destination[seq];
                new_proteins.Add(
                    new Protein(
                        seq,
                        destination_protein.Accession,
                        gene_names : source_protein.GeneNames.ToList(),
                        oneBasedModifications : source_protein.OneBasedPossibleLocalizedModifications, 
                        proteolysisProducts : source_protein.ProteolysisProducts.ToList(),
                        name : destination_protein.Name, 
                        full_name : destination_protein.FullName,
                        isDecoy : destination_protein.IsDecoy,
                        isContaminant : destination_protein.IsContaminant, 
                        databaseReferences : source_protein.DatabaseReferences.ToList(),
                        sequenceVariations : source_protein.SequenceVariations.ToList()));
            }

            List<string> destination_only = dict_destination.Keys.Except(dict_source.Keys).ToList();
            List<string> source_only = dict_source.Keys.Except(dict_destination.Keys).ToList();
            return new_proteins;
        }
    }
}
