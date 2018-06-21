using Bio.VCF;
using Proteogenomics;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

namespace WorkflowLayer
{
    public static class ProteinAnnotation
    {
        /// <summary>
        /// Gets UniProt ptmlist
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public static List<ModificationWithLocation> GetUniProtMods(string spritzDirectory)
        {
            Loaders.LoadElements(Path.Combine(spritzDirectory, "elements.dat"));
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(spritzDirectory, "PSI-MOD.obo.xml"));
            return Loaders.LoadUniprot(Path.Combine(spritzDirectory, "ptmlist.txt"), Loaders.GetFormalChargesDictionary(psiModDeserialized)).ToList();
        }

        /// <summary>
        /// Takes in a protein XML and completes variant annotations by sorting them and annotating genotypes
        /// </summary>
        /// <param name="proteinXmlPath"></param>
        /// <param name="vcfPath"></param>
        /// <returns>Path of newly written protein XML</returns>
        public static string CompleteVariantAnnotations(string spritzDirectory, string proteinXmlPath, string vcfPath, Genome genome)
        {
            List<Protein> outprots = new List<Protein>();
            var proteins = ProteinDbLoader.LoadProteinXML(proteinXmlPath, true, DecoyType.None, GetUniProtMods(spritzDirectory), false, null, out var un);
            var variants = new VCFParser(vcfPath).Select(v => new Variant(null, v, genome));
            IntervalForest forest = new IntervalForest(variants);
            foreach (var p in proteins)
            {
                Protein outprot = p;
                if (p.SequenceVariations.Count() > 0)
                {
                    List<SequenceVariation> completedVariations = p.SequenceVariations.Select(v => CompleteVariationAnnotation(v, forest))
                        .Where(v => v != null)
                        .OrderBy(v => v.OneBasedBeginPosition).ToList();
                    outprot = new Protein(
                        p.BaseSequence,
                        p.Accession,
                        organism: p.Organism,
                        name: p.Name,
                        full_name: p.FullName,
                        isDecoy: p.IsDecoy,
                        isContaminant: p.IsContaminant,
                        sequenceVariations: completedVariations,

                        // combine these
                        gene_names: p.GeneNames.ToList(),

                        // transfer these
                        oneBasedModifications: p.OneBasedPossibleLocalizedModifications,
                        proteolysisProducts: p.ProteolysisProducts.ToList(),
                        databaseReferences: p.DatabaseReferences.ToList(),
                        disulfideBonds: p.DisulfideBonds.ToList());
                }
                outprots.Add(outprot);
            }
            string outProteinXmlPath = Path.Combine(Path.GetDirectoryName(proteinXmlPath), Path.GetFileNameWithoutExtension(proteinXmlPath) + ".spritzAnn.xml");
            ProteinDbWriter.WriteXmlDatabase(null, outprots, outProteinXmlPath);
            return outProteinXmlPath;
        }

        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="proteogenomicProteins"></param>
        /// <param name="uniprotProteins"></param>
        /// <returns></returns>
        public static List<Protein> CombineAndAnnotateProteins(List<Protein> uniprotProteins, List<Protein> proteogenomicProteins)
        {
            List<Protein> newProteins = new List<Protein>();
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

                if (uniprot != null && uniprot.BaseSequence != pgProtein.BaseSequence) { throw new ArgumentException("not all proteins have the same sequence"); }

                // keep all sequence variations
                newProteins.Add(new Protein(
                    pgProtein.BaseSequence,
                    uniprot != null ? uniprot.Accession : pgProtein.Accession, //comma-separated
                    organism: uniprot != null ? uniprot.Organism : pgProtein.Organism,
                    name: pgProtein.Name, //comma-separated
                    full_name: uniprot != null ? uniprot.FullName : pgProtein.FullName, //comma-separated
                    isDecoy: pgProtein.IsDecoy,
                    isContaminant: pgProtein.IsContaminant,
                    sequenceVariations: pgProtein.SequenceVariations.OrderBy(v => v.OneBasedBeginPosition).ToList(),

                    // combine these
                    gene_names: (uniprot != null ? uniprot.GeneNames : new List<Tuple<string, string>>()).Concat(pgProtein.GeneNames).ToList(),

                    // transfer these
                    oneBasedModifications: uniprot != null ? uniprot.OneBasedPossibleLocalizedModifications : new Dictionary<int, List<Modification>>(),
                    proteolysisProducts: uniprot != null ? uniprot.ProteolysisProducts.ToList() : new List<ProteolysisProduct>(),
                    databaseReferences: uniprot != null ? uniprot.DatabaseReferences.ToList() : new List<DatabaseReference>(),
                    disulfideBonds: uniprot != null ? uniprot.DisulfideBonds.ToList() : new List<DisulfideBond>()
                ));
            }

            return newProteins;
        }

        /// <summary>
        /// Completes a protein sequence variant annotation with variants from a VCF file
        /// </summary>
        /// <param name="variation"></param>
        /// <param name="forestWithVariants"></param>
        /// <returns></returns>
        private static SequenceVariation CompleteVariationAnnotation(SequenceVariation variation, IntervalForest forestWithVariants)
        {
            string[] fields = variation.Description.Split(new[] { ' ', '\t' });
            if (fields.Length < 5) // this occurs for most uniprot variations, which should be ignored in this use case
            {
                return null;
            }
            string chromosomeName = fields[0];
            if (!int.TryParse(fields[1], out int start) || !forestWithVariants.Forest.TryGetValue(chromosomeName, out IntervalTree tree))
            {
                return null;
            }

            List<string[]> variantLines = tree.Stab(start).Distinct().Select(i => (i as Variant).VariantContext.VcfLine.Split('\t'))
                .Where(line => int.TryParse(line[1], out int vcfLineStart) && start == vcfLineStart)
                .ToList();
            string referenceAlleleString = fields[3];
            string alternateAlleleString = fields[4];
            string[] closestVariant = variantLines
                .Where(line => line[3] == referenceAlleleString && line[4] == alternateAlleleString
                    || alternateAlleleString == "" && line[3].EndsWith(referenceAlleleString) // deletion in SnpEff simplified format
                    || referenceAlleleString == "" && alternateAlleleString.EndsWith(line[4])) // insertion in SnpEff simplified format
                .FirstOrDefault();
            variation.Description = String.Join("\t", closestVariant);
            return variation;
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

            return new Protein(
                seq.First(),
                String.Join(",", proteinsWithSameSequence.Select(p => p.Accession)),
                organism: proteinsWithSameSequence.First().Organism,
                name: String.Join(",", proteinsWithSameSequence.Select(p => p.Name)),
                full_name: String.Join(",", proteinsWithSameSequence.Select(p => p.FullName)),
                isDecoy: proteinsWithSameSequence.All(p => p.IsDecoy),
                isContaminant: proteinsWithSameSequence.Any(p => p.IsContaminant),
                sequenceVariations: proteinsWithSameSequence.SelectMany(p => p.SequenceVariations).ToList(),
                gene_names: proteinsWithSameSequence.SelectMany(p => p.GeneNames).ToList(),
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
            Dictionary<int, HashSet<Modification>> result = new Dictionary<int, HashSet<Modification>>();
            if (proteins != null)
            {
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