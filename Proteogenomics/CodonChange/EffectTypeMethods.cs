using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public static class EffectTypeMethods
    {
        public static readonly Dictionary<EffectType, EffectImpact> EffectDictionary = new Dictionary<EffectType, EffectImpact>
        {
            // High impact
            // Order: Highest impact first
            { EffectType.CHROMOSOME_LARGE_DELETION, EffectImpact.HIGH }, //
            { EffectType.CHROMOSOME_LARGE_INVERSION, EffectImpact.HIGH }, //
            { EffectType.CHROMOSOME_LARGE_DUPLICATION, EffectImpact.HIGH }, //
            { EffectType.GENE_REARRANGEMENT, EffectImpact.HIGH }, //
            { EffectType.GENE_DELETED, EffectImpact.HIGH }, //
            { EffectType.TRANSCRIPT_DELETED, EffectImpact.HIGH }, //
            { EffectType.EXON_DELETED, EffectImpact.HIGH }, //
            { EffectType.EXON_DELETED_PARTIAL, EffectImpact.HIGH }, //
            { EffectType.GENE_FUSION, EffectImpact.HIGH }, //
            { EffectType.GENE_FUSION_REVERESE, EffectImpact.HIGH }, //
            { EffectType.GENE_FUSION_HALF, EffectImpact.HIGH }, //
            { EffectType.FRAME_SHIFT, EffectImpact.HIGH }, //
            { EffectType.STOP_GAINED, EffectImpact.HIGH }, //
            { EffectType.STOP_LOST, EffectImpact.HIGH }, //
            { EffectType.START_LOST, EffectImpact.HIGH }, //
            { EffectType.SPLICE_SITE_ACCEPTOR, EffectImpact.HIGH }, //
            { EffectType.SPLICE_SITE_DONOR, EffectImpact.HIGH }, //
            { EffectType.RARE_AMINO_ACID, EffectImpact.HIGH }, //
            { EffectType.EXON_DUPLICATION, EffectImpact.HIGH }, //
            { EffectType.EXON_DUPLICATION_PARTIAL, EffectImpact.HIGH }, //
            { EffectType.EXON_INVERSION, EffectImpact.HIGH }, //
            { EffectType.EXON_INVERSION_PARTIAL, EffectImpact.HIGH }, //
            { EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS, EffectImpact.HIGH }, //
            { EffectType.PROTEIN_STRUCTURAL_INTERACTION_LOCUS, EffectImpact.HIGH }, //

            // Moderate impact
            // Order: Highest impact first
            // Note: Method Codon.effect() relies on this order for effect
            //       replacement (when 'allowReplace = true')
            { EffectType.NON_SYNONYMOUS_CODING, EffectImpact.MODERATE }, //
            { EffectType.GENE_DUPLICATION, EffectImpact.MODERATE }, //
            { EffectType.TRANSCRIPT_DUPLICATION, EffectImpact.MODERATE }, //
            { EffectType.UTR_5_DELETED, EffectImpact.MODERATE }, //
            { EffectType.UTR_3_DELETED, EffectImpact.MODERATE }, //
            { EffectType.SPLICE_SITE_BRANCH_U12, EffectImpact.MODERATE }, //
            { EffectType.GENE_INVERSION, EffectImpact.MODERATE }, //
            { EffectType.TRANSCRIPT_INVERSION, EffectImpact.MODERATE }, //
            { EffectType.CODON_INSERTION, EffectImpact.MODERATE }, //
            { EffectType.CODON_CHANGE_PLUS_CODON_INSERTION, EffectImpact.MODERATE }, //
            { EffectType.CODON_DELETION, EffectImpact.MODERATE }, //
            { EffectType.CODON_CHANGE_PLUS_CODON_DELETION, EffectImpact.MODERATE }, //

            // Low impact
            // Order: Highest impact first
            { EffectType.NON_SYNONYMOUS_STOP, EffectImpact.LOW }, //
            { EffectType.NON_SYNONYMOUS_START, EffectImpact.LOW }, //
            { EffectType.SPLICE_SITE_REGION, EffectImpact.LOW }, //
            { EffectType.SPLICE_SITE_BRANCH, EffectImpact.LOW }, //
            { EffectType.SYNONYMOUS_CODING, EffectImpact.LOW }, //
            { EffectType.SYNONYMOUS_START, EffectImpact.LOW }, //
            { EffectType.SYNONYMOUS_STOP, EffectImpact.LOW }, //
            { EffectType.CODON_CHANGE, EffectImpact.LOW }, //
            { EffectType.START_GAINED, EffectImpact.LOW }, //
            { EffectType.MOTIF, EffectImpact.LOW }, //
            { EffectType.MOTIF_DELETED, EffectImpact.LOW }, //
            { EffectType.FEATURE_FUSION, EffectImpact.LOW }, //

            // Modifiers
            // Order: Highest impact first
            { EffectType.UTR_5_PRIME, EffectImpact.MODIFIER }, //
            { EffectType.UTR_3_PRIME, EffectImpact.MODIFIER }, //
            { EffectType.REGULATION, EffectImpact.MODIFIER }, //
            { EffectType.MICRO_RNA, EffectImpact.MODIFIER }, //
            { EffectType.UPSTREAM, EffectImpact.MODIFIER }, //
            { EffectType.DOWNSTREAM, EffectImpact.MODIFIER }, //
            { EffectType.NEXT_PROT, EffectImpact.MODIFIER }, //
            { EffectType.INTRON_CONSERVED, EffectImpact.MODIFIER }, //
            { EffectType.INTRON, EffectImpact.MODIFIER }, //
            { EffectType.INTRAGENIC, EffectImpact.MODIFIER }, //
            { EffectType.INTERGENIC_CONSERVED, EffectImpact.MODIFIER }, //
            { EffectType.INTERGENIC, EffectImpact.MODIFIER }, //
            { EffectType.CDS, EffectImpact.MODIFIER }, //
            { EffectType.EXON, EffectImpact.MODIFIER }, //
            { EffectType.TRANSCRIPT, EffectImpact.MODIFIER }, //
            { EffectType.GENE, EffectImpact.MODIFIER }, //
            { EffectType.SEQUENCE, EffectImpact.MODIFIER }, //
            { EffectType.CHROMOSOME_ELONGATION, EffectImpact.MODIFIER }, //
            { EffectType.CUSTOM, EffectImpact.MODIFIER }, //
            { EffectType.CHROMOSOME, EffectImpact.MODIFIER }, //
            { EffectType.GENOME, EffectImpact.MODIFIER }, //
            { EffectType.NONE, EffectImpact.MODIFIER } //
        };

        public static string GetEffectTypeString(bool useFirstEffect, List<EffectType> EffectTypes)
        {
            if (EffectTypes == null) { return ""; }

            // Show all effects
            StringBuilder sb = new StringBuilder();
            EffectTypes.Sort();

            // More than one effect? Check for repeats
            HashSet<string> added = EffectTypes.Count > 1 && !useFirstEffect ? new HashSet<string>() : null;

            // Create string
            foreach (EffectType et in EffectTypes)
            {
                string eff = et.ToString();

                // Make sure we don't add the same effect twice
                if (added == null || added.Add(eff))
                {
                    sb.Append(sb.Length > 0 ? " " : eff);
                }

                // Only use first effect?
                if (useFirstEffect)
                {
                    return sb.ToString();
                }
            }

            return sb.ToString();
        }

        public static bool HasEffectType(EffectType effectType, IEnumerable<EffectType> effectTypes)
        {
            return effectTypes.Any(t => effectType == t);
        }

        public static bool IsExon(Interval marker, IEnumerable<EffectType> effectTypes)
        {
            return marker is Exon || HasEffectType(EffectType.EXON_DELETED, effectTypes);
        }

        public static bool IsIntergenic(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.INTERGENIC, effectTypes) || HasEffectType(EffectType.INTERGENIC_CONSERVED, effectTypes);
        }

        public static bool IsIntron(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.INTRON, effectTypes) || HasEffectType(EffectType.INTRON_CONSERVED, effectTypes);
        }

        public static bool IsMotif(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.MOTIF, effectTypes);
        }

        public static bool IsNextProt(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.NEXT_PROT, effectTypes);
        }

        public static bool IsRegulation(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.REGULATION, effectTypes);
        }

        public static bool IsSpliceSite(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.SPLICE_SITE_DONOR, effectTypes) //
                    || HasEffectType(EffectType.SPLICE_SITE_ACCEPTOR, effectTypes) //
                    || HasEffectType(EffectType.SPLICE_SITE_REGION, effectTypes) //
                    || HasEffectType(EffectType.SPLICE_SITE_BRANCH, effectTypes) //
                    || HasEffectType(EffectType.SPLICE_SITE_BRANCH_U12, effectTypes); //
        }

        public static bool IsSpliceSiteCore(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.SPLICE_SITE_DONOR, effectTypes) //
                    || HasEffectType(EffectType.SPLICE_SITE_ACCEPTOR, effectTypes); //
        }

        public static bool IsSpliceSiteRegion(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.SPLICE_SITE_REGION, effectTypes);
        }

        public static bool IsUtr3(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.UTR_3_PRIME, effectTypes) || HasEffectType(EffectType.UTR_3_DELETED, effectTypes);
        }

        public static bool IsUtr5(IEnumerable<EffectType> effectTypes)
        {
            return HasEffectType(EffectType.UTR_5_PRIME, effectTypes) || HasEffectType(EffectType.UTR_5_DELETED, effectTypes);
        }

        public static EffectType GetGeneRegion(EffectType type)
        {
            switch (type)
            {
                case EffectType.NONE:
                case EffectType.CHROMOSOME:
                case EffectType.CHROMOSOME_LARGE_DELETION:
                case EffectType.CHROMOSOME_LARGE_DUPLICATION:
                case EffectType.CHROMOSOME_LARGE_INVERSION:
                case EffectType.CHROMOSOME_ELONGATION:
                case EffectType.CUSTOM:
                case EffectType.SEQUENCE:
                    return EffectType.CHROMOSOME;

                case EffectType.INTERGENIC:
                case EffectType.INTERGENIC_CONSERVED:
                case EffectType.FEATURE_FUSION:
                    return EffectType.INTERGENIC;

                case EffectType.UPSTREAM:
                    return EffectType.UPSTREAM;

                case EffectType.UTR_5_PRIME:
                case EffectType.UTR_5_DELETED:
                case EffectType.START_GAINED:
                    return EffectType.UTR_5_PRIME;

                case EffectType.SPLICE_SITE_ACCEPTOR:
                    return EffectType.SPLICE_SITE_ACCEPTOR;

                case EffectType.SPLICE_SITE_BRANCH_U12:
                case EffectType.SPLICE_SITE_BRANCH:
                    return EffectType.SPLICE_SITE_BRANCH;

                case EffectType.SPLICE_SITE_DONOR:
                    return EffectType.SPLICE_SITE_DONOR;

                case EffectType.SPLICE_SITE_REGION:
                    return EffectType.SPLICE_SITE_REGION;

                case EffectType.TRANSCRIPT_DELETED:
                case EffectType.TRANSCRIPT_DUPLICATION:
                case EffectType.TRANSCRIPT_INVERSION:
                case EffectType.INTRAGENIC:
                case EffectType.NEXT_PROT:
                case EffectType.TRANSCRIPT:
                case EffectType.CDS:
                    return EffectType.TRANSCRIPT;

                case EffectType.GENE:
                case EffectType.GENE_DELETED:
                case EffectType.GENE_DUPLICATION:
                case EffectType.GENE_FUSION:
                case EffectType.GENE_FUSION_HALF:
                case EffectType.GENE_FUSION_REVERESE:
                case EffectType.GENE_INVERSION:
                case EffectType.GENE_REARRANGEMENT:
                    return EffectType.GENE;

                case EffectType.EXON:
                case EffectType.EXON_DELETED:
                case EffectType.EXON_DELETED_PARTIAL:
                case EffectType.EXON_DUPLICATION:
                case EffectType.EXON_DUPLICATION_PARTIAL:
                case EffectType.EXON_INVERSION:
                case EffectType.EXON_INVERSION_PARTIAL:
                case EffectType.NON_SYNONYMOUS_START:
                case EffectType.NON_SYNONYMOUS_CODING:
                case EffectType.SYNONYMOUS_CODING:
                case EffectType.SYNONYMOUS_START:
                case EffectType.FRAME_SHIFT:
                case EffectType.CODON_CHANGE:
                case EffectType.CODON_INSERTION:
                case EffectType.CODON_CHANGE_PLUS_CODON_INSERTION:
                case EffectType.CODON_DELETION:
                case EffectType.CODON_CHANGE_PLUS_CODON_DELETION:
                case EffectType.START_LOST:
                case EffectType.STOP_GAINED:
                case EffectType.SYNONYMOUS_STOP:
                case EffectType.NON_SYNONYMOUS_STOP:
                case EffectType.STOP_LOST:
                case EffectType.RARE_AMINO_ACID:
                case EffectType.PROTEIN_PROTEIN_INTERACTION_LOCUS:
                case EffectType.PROTEIN_STRUCTURAL_INTERACTION_LOCUS:
                    return EffectType.EXON;

                case EffectType.INTRON:
                case EffectType.INTRON_CONSERVED:
                    return EffectType.INTRON;

                case EffectType.UTR_3_PRIME:
                case EffectType.UTR_3_DELETED:
                    return EffectType.UTR_3_PRIME;

                case EffectType.DOWNSTREAM:
                    return EffectType.DOWNSTREAM;

                case EffectType.REGULATION:
                    return EffectType.REGULATION;

                case EffectType.MOTIF:
                case EffectType.MOTIF_DELETED:
                    return EffectType.MOTIF;

                case EffectType.MICRO_RNA:
                    return EffectType.MICRO_RNA;

                case EffectType.GENOME:
                    return EffectType.GENOME;

                default:
                    throw new ArgumentException("Unknown gene region for effect type: '" + type + "'");
            }
        }
    }
}