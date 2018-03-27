using Bio;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace Proteogenomics
{
    public class Gene :
        Interval
    {
        private MetadataListItem<List<string>> Metadata { get; set; }

        public string ID { get; set; }
        public Chromosome Chromosome { get; set; }
        public List<Transcript> Transcripts { get; set; } = new List<Transcript>();
        public IntervalTree TranscriptTree { get; set; }

        public Gene(string ID, Chromosome chromosome, string strand, long oneBasedStart, long oneBasedEnd, MetadataListItem<List<string>> metadata)
            : base(chromosome, chromosome.Sequence.ID, strand, oneBasedStart, oneBasedEnd)
        {
            this.ID = ID;
            this.Chromosome = chromosome;
            this.Metadata = metadata;
        }

        public List<Protein> Translate(bool translateCodingDomains, bool includeVariants, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            List<Protein> proteins = new List<Protein>();
            foreach (Transcript t in Transcripts)
            {
                if (incompleteTranscriptAccessions.Contains(t.ProteinID)) { continue; }
                lock (proteins) { proteins.Add(t.Protein(selenocysteineContaining)); }
            }
            return proteins;
        }

        public List<Protein> TranslateUsingAnnotatedStartCodons(GeneModel genesWithCodingDomainSequences, Dictionary<string, string> selenocysteineContaining, bool includeVariants, int min_length)
        {
            int indexBinSize = 100000;
            Dictionary<Tuple<string, string, long>, List<CDS>> binnedCodingStarts = new Dictionary<Tuple<string, string, long>, List<CDS>>();
            foreach (CDS x in genesWithCodingDomainSequences.Genes.SelectMany(g => g.Transcripts).Select(t => t.CodingDomainSequences.FirstOrDefault()).Where(x => x != null))
            {
                Tuple<string, string, long> key = new Tuple<string, string, long>(x.ChromosomeID, x.Strand, x.OneBasedStart / indexBinSize * indexBinSize);
                binnedCodingStarts.TryGetValue(key, out List<CDS> xs);
                if (xs == null) { binnedCodingStarts.Add(key, new List<CDS> { x }); }
                else { binnedCodingStarts[key].Add(x); }
            }

            List<Protein> proteins = new List<Protein>();
            Parallel.ForEach(Transcripts, t =>
            {
                IEnumerable<Protein> p = t.TranslateUsingAnnotatedStartCodons(binnedCodingStarts, selenocysteineContaining, indexBinSize, min_length, includeVariants);
                if (p != null)
                {
                    lock (proteins) proteins.AddRange(p);
                }
            });
            return proteins;
        }
    }
}