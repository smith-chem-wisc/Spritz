using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Bio;
using Proteomics;

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
                Protein p = t.translate(translateCDS, includeVariants);
                if (p != null)
                    lock (proteins) proteins.Add(p);
            });
            return proteins;
        }

        public List<Protein> translate(GeneModel genesWithCDS, int min_length, bool includeVariants)
        {
            List<Protein> proteins = new List<Protein>();
            Parallel.ForEach(transcripts, t =>
            {
                Protein p = t.translate(genesWithCDS, min_length, includeVariants);
                if (p != null)
                    lock (proteins) proteins.Add(p);
            });
            return proteins;
        }

        #endregion Public Methods

    }
}
