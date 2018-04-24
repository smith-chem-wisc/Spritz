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
        /// <summary>
        /// Constructing a gene using Bio object
        /// </summary>
        /// <param name="id"></param>
        /// <param name="chromosome"></param>
        /// <param name="strand"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="metadata"></param>
        public Gene(string id, Chromosome chromosome, string strand, long oneBasedStart, long oneBasedEnd)
            : base(chromosome, chromosome.Sequence.ID, strand, oneBasedStart, oneBasedEnd, null)
        {
            ID = id;
            Chromosome = chromosome;
        }

        public string ID { get; set; }
        public Chromosome Chromosome { get; set; }
        public List<Transcript> Transcripts { get; set; } = new List<Transcript>();
        public IntervalTree TranscriptTree { get; set; } = new IntervalTree();

        public List<Protein> Translate(bool translateCodingDomains, HashSet<string> incompleteTranscriptAccessions = null, Dictionary<string, string> selenocysteineContaining = null)
        {
            incompleteTranscriptAccessions = incompleteTranscriptAccessions ?? new HashSet<string>();
            List<Protein> proteins = new List<Protein>();
            foreach (Transcript t in translateCodingDomains ? Transcripts.Where(t => t.IsProteinCoding()) : Transcripts)
            {
                if (incompleteTranscriptAccessions.Contains(t.ProteinID)) { continue; }
                lock (proteins) { proteins.Add(t.Protein(selenocysteineContaining)); }
            }
            return proteins;
        }

        public List<Protein> TranslateUsingAnnotatedStartCodons(GeneModel genesWithCodingDomainSequences, Dictionary<string, string> selenocysteineContaining, int min_length)
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
                IEnumerable<Protein> p = t.TranslateUsingAnnotatedStartCodons(binnedCodingStarts, selenocysteineContaining, indexBinSize, min_length);
                if (p != null)
                {
                    lock (proteins) proteins.AddRange(p);
                }
            });
            return proteins;
        }
    }
}