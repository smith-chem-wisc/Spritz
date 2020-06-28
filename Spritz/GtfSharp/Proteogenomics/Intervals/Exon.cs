using Bio;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Exon :
        IntervalSequence
    {
        /// <summary>
        /// Construct an exon
        /// </summary>
        /// <param name="parent"></param>
        /// <param name="Sequence"></param>
        /// <param name="oneBasedStart"></param>
        /// <param name="oneBasedEnd"></param>
        /// <param name="chromID"></param>
        /// <param name="strand"></param>
        public Exon(Transcript parent, ISequence Sequence, string source, long oneBasedStart, long oneBasedEnd, string chromID, string strand, MetadataListItem<List<string>> featureMetadata)
            : base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd, Sequence)
        {
            FeatureMetadata = featureMetadata;
        }

        /// <summary>
        /// Copy an exon
        /// </summary>
        /// <param name="x"></param>
        public Exon(Exon x)
            : this(x.Parent as Transcript, x.Sequence, x.Source, x.OneBasedStart, x.OneBasedEnd, x.ChromosomeID, x.Strand, x.FeatureMetadata)
        {
        }

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "exon";

        public MetadataListItem<List<string>> FeatureMetadata { get; private set; }

        public override string GetGtfAttributes()
        {
            var attributes = GeneModel.SplitAttributes(FeatureMetadata.FreeText);
            List<Tuple<string, string>> attributeSubsections = new List<Tuple<string, string>>();

            string exonIdLabel = "exon_id";
            bool hasExonId = attributes.TryGetValue(exonIdLabel, out string exonId);
            if (hasExonId) { attributeSubsections.Add(new Tuple<string, string>(exonIdLabel, exonId)); }

            string exonVersionLabel = "exon_version";
            bool hasExonVersion = attributes.TryGetValue(exonVersionLabel, out string exonVersion);
            if (hasExonVersion) { attributeSubsections.Add(new Tuple<string, string>(exonVersionLabel, exonVersion)); }

            string exonNumberLabel = "exon_number";
            string exonNumber = (Parent as Transcript).Exons.Count(x => x.OneBasedStart <= OneBasedStart).ToString();
            attributeSubsections.Add(new Tuple<string, string>(exonNumberLabel, exonNumber));

            return Parent.GetGtfAttributes() + " " + String.Join(" ", attributeSubsections.Select(x => x.Item1 + " \"" + x.Item2 + "\";"));
        }
    }
}