using Bio;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Linq;

namespace Genomics
{
    public class Gene
    {

        #region Private Properties

        private MetadataListItem<List<string>> metadata { get; set; }

        #endregion

        #region Public Properties

        public string ID { get; set; }
        public string ChromID { get; set; }
        public List<Transcript> transcripts { get; set; } = new List<Transcript>();
        public HashSet<Exon> exons { get; set; } = new HashSet<Exon>();
        public HashSet<Exon> CDS { get; set; } = new HashSet<Exon>();

        #endregion

        #region Public Constructors

        public Gene(string ID, string ChromID, MetadataListItem<List<string>> metadata)
        {
            this.ID = ID;
            this.ChromID = ChromID;
            this.metadata = metadata;
        }

        #endregion Public Constructors

        #region Public Methods

        public List<Protein> translate(bool translateCDS, bool includeVariants)
        {
            List<Protein> proteins = new List<Protein>();
            Parallel.ForEach(transcripts, t =>
            {
                IEnumerable<Protein> p = t.translate(translateCDS, includeVariants);
                if (p != null)
                    lock (proteins) proteins.AddRange(p);
            });
            return proteins;
        }

        public List<Protein> translate(GeneModel genesWithCDS, int min_length, bool includeVariants)
        {
            int bin_size = 100000;
            Dictionary<Tuple<string, string, long>, List<Exon>> chr_index_CDS = new Dictionary<Tuple<string, string, long>, List<Exon>>();
            foreach (Exon x in genesWithCDS.genes.SelectMany(g => g.transcripts).SelectMany(t => t.CDS))
            {
                Tuple<string, string, long> key = new Tuple<string, string, long>(x.ChromID, x.Strand, x.OneBasedStart / bin_size * bin_size);
                chr_index_CDS.TryGetValue(key, out List<Exon> xs);
                if (xs == null) chr_index_CDS.Add(key, new List<Exon> { x });
                else chr_index_CDS[key].Add(x);
            }

            List<Protein> proteins = new List<Protein>();
            Parallel.ForEach(transcripts, t =>
            {
                IEnumerable<Protein> p = t.translate(chr_index_CDS, bin_size, min_length, includeVariants);
                if (p != null)
                    lock (proteins) proteins.AddRange(p);
            });
            return proteins;
        }

        #endregion Public Methods

    }
}
