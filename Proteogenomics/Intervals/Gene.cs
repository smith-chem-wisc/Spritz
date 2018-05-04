using Proteomics;
using System.Collections.Generic;
using System.Linq;

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
    }
}