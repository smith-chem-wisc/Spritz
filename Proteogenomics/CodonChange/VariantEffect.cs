using Bio;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class VariantEffect
        : IComparable<VariantEffect>
    {
        // Don't show codon change sequences that are too long
        public static readonly int MAX_CODON_SEQUENCE_LEN = 100;

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

        public VariantEffect(Variant variant)
        {
            Variant = variant;
            EffectTypes = new List<EffectType>();
            EffectImpacts = new List<EffectImpact>();
        }

        public VariantEffect(Variant variant, Interval marker, EffectType effectType, EffectImpact effectImpact, string codonsOld, string codonsNew, int codonNum, int codonIndex, long cDnaPos)
        {
            Variant = variant;
            EffectTypes = new List<EffectType>();
            EffectImpacts = new List<EffectImpact>();
            Set(marker, effectType, effectImpact, "");
            SetCodons(codonsOld, codonsNew, codonNum, codonIndex);
            CDnaPos = cDnaPos;
        }

        public Variant Variant { get; private set; }

        public List<EffectType> EffectTypes { get; private set; }

        public string CodonsRef { get; private set; } = ""; // Codon change information

        public string CodonsAlt { get; private set; } = ""; // Codon change information

        public long CDnaPos { get; private set; } = -1; // Position in cDNA

        public int CodonNum { get; private set; } = -1; // Codon number (negative number mens 'information not available')

        public int CodonIndex { get; private set; } = -1; // Index within a codon (negative number mens 'information not available')

        public int CodonDegeneracy { get; private set; } = -1; // Codon degeneracy (negative number mens 'information not available')

        public string ReferenceAA { get; set; } = "";

        public string AlternateAA { get; private set; } = ""; // Amino acid changes

        protected List<EffectImpact> EffectImpacts { get; set; }

        protected Interval Marker { get; set; }

        protected HashSet<ErrorWarningType> Error { get; set; } = new HashSet<ErrorWarningType>();

        protected HashSet<ErrorWarningType> Warning { get; set; } = new HashSet<ErrorWarningType>();

        protected string Message { get; set; } = ""; // Any message, warning or error?

        protected string CodonsAroundOld { get; private set; } = ""; // Codons around

        protected string CodonsAroundNew { get; private set; } = "";

        protected long Distance = -1; // Distance metric

        protected string AroundOldAAs { get; set; } = "";

        protected string AroundNewAAs { get; set; } = ""; // Amino acids around

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

        public void AddEffect(EffectType effectType)
        {
            AddEffectType(effectType);
            AddEffectImpact(EffectDictionary[effectType]);
        }

        public void AddEffectImpact(EffectImpact effectImpact)
        {
            EffectImpacts.Add(effectImpact);
        }

        public void AddEffectType(EffectType effectType)
        {
            EffectTypes.Add(effectType);
        }

        /// <summary>
        /// Create a string for codon effect
        /// @param showAaChange : If true, include codon change, biotype, etc.
        /// </summary>
        /// <param name="showAaChange"></param>
        /// <param name="showBioType"></param>
        /// <param name="useSeqOntology"></param>
        /// <param name="useFirstEffect"></param>
        /// <returns></returns>
        private string CodonEffect(bool showAaChange, bool showBioType, bool useFirstEffect)
        {
            if ((Marker == null) || (CodonNum < 0)) { return ""; }

            if (!showAaChange) { return GetEffectTypeString(useFirstEffect); }

            StringBuilder sb = new StringBuilder();
            sb.Append(GetEffectTypeString(useFirstEffect));
            sb.Append("(");
            sb.Append(GetAaChange());
            sb.Append(")");

            return sb.ToString();
        }

        /// <summary>
        /// Compare to another variantEffect
        /// </summary>
        /// <param name="variantEffect"></param>
        /// <returns></returns>
        public int CompareTo(VariantEffect varEffOther)
        {
            // Sort by impact
            int comp = GetEffectImpact().CompareTo(varEffOther.GetEffectImpact());
            if (comp != 0) return comp;

            // Sort by effect
            comp = GetEffectType().CompareTo(varEffOther.GetEffectType());
            if (comp != 0) return comp;

            //---
            // Transcript based comparisons
            //---
            Transcript trThis = GetTranscript();
            Transcript trOther = varEffOther.GetTranscript();

            // This transcript data
            //int tslThis = TranscriptSupportLevel.TSL_NULL_VALUE;
            //int canonThis = 0;
            long startThis = 0;
            long endThis = 0;
            string idThis = null;
            string chrThis = null;
            if (trThis != null)
            {
                //tslThis = TranscriptSupportLevel.tsl(trThis.getTranscriptSupportLevel());
                //canonThis = trThis.isCanonical() ? 1 : 0;
                idThis = trThis.ID;
                chrThis = trThis.ChromosomeID;
                startThis = trThis.OneBasedStart;
                endThis = trThis.OneBasedEnd;
            }

            // Other transcript data
            //int tslOther = TranscriptSupportLevel.TSL_NULL_VALUE;
            //int canonOther = 0;
            long startOther = 0;
            long endOther = 0;
            string idOther = null;
            string chrOther = null;
            if (trOther != null)
            {
                //tslOther = TranscriptSupportLevel.tsl(trOther.getTranscriptSupportLevel());
                //canonOther = trOther.isCanonical() ? 1 : 0;
                idOther = trOther.ID;
                chrOther = trOther.ChromosomeID;
                startOther = trOther.OneBasedStart;
                endOther = trOther.OneBasedEnd;
            }

            // Compare by TSL: Smaller first
            //comp = tslThis - tslOther;
            //if (comp != 0) return comp;

            // Compare by canonical transcript: Larger first
            //comp = canonOther - canonThis;
            //if (comp != 0) return comp;

            // Compare genomic coordinate
            comp = CompareNull(chrThis, chrOther);
            if (comp != 0) { return comp; }

            comp = startThis.CompareTo(startOther);
            if (comp != 0) { return comp; }

            comp = endThis.CompareTo(endOther);
            if (comp != 0) { return comp; }

            // Compare IDs
            comp = CompareNull(idThis, idOther);
            if (comp != 0) { return comp; }

            //---
            // Marker based comparisons
            //---
            Interval mThis = GetMarker();
            Interval mOther = varEffOther.GetMarker();

            startThis = 0;
            endThis = 0;
            idThis = null;
            chrThis = null;
            if (mThis != null)
            {
                //idThis = mThis.getId();
                chrThis = mThis.ChromosomeID;
                startThis = mThis.OneBasedStart;
                endThis = mThis.OneBasedEnd;
            }

            startOther = 0;
            endOther = 0;
            idOther = null;
            chrOther = null;
            if (mOther != null)
            {
                //idOther = mOther.getId();
                chrOther = mOther.ChromosomeID;
                startOther = mOther.OneBasedStart;
                endOther = mOther.OneBasedEnd;
            }

            // Compare genomic coordinate
            comp = CompareNull(chrThis, chrOther);
            if (comp != 0) { return comp; }

            comp = startThis.CompareTo(startOther);
            if (comp != 0) { return comp; }

            comp = endThis.CompareTo(endOther);
            if (comp != 0) { return comp; }

            // Compare IDs
            comp = CompareNull(idThis, idOther);
            if (comp != 0) { return comp; }

            //---
            // Variant based comparison
            //---
            return Variant.CompareTo(varEffOther.Variant);
        }

        /// <summary>
        /// Show a string with overall effect
        /// </summary>
        /// <param name="shortFormat"></param>
        /// <param name="showAaChange"></param>
        /// <param name="showBioType"></param>
        /// <param name="useSeqOntology"></param>
        /// <param name="useFirstEffect"></param>
        /// <returns></returns>
        public string Effect(bool shortFormat, bool showAaChange, bool showBioType, bool useSeqOntology, bool useFirstEffect)
        {
            string e = "";
            string codonEffect = CodonEffect(showAaChange, showBioType, useFirstEffect); // Codon effect

            // Create effect string
            if (codonEffect != "")
            {
                e = codonEffect;
            }
            else if (IsIntergenic() || IsIntron() || IsSpliceSite())
            {
                e = GetEffectTypeString(useFirstEffect);
            }
            else if (Message != "")
            {
                e = GetEffectTypeString(useFirstEffect) + ": " + Message;
            }
            else if (Marker == null)
            {
                // There are cases when no marker is associated (e.g. "Out of chromosome", "No such chromosome", etc.)
                e = GetEffectTypeString(useFirstEffect);
            }
            else
            {
                e = GetEffectTypeString(useFirstEffect); // + ": " + marker.getId();
            }

            if (shortFormat)
            {
                e = e.Split(':')[0];
            }

            return e;
        }

        /// <summary>
        /// Amino acid change string (HGVS style)
        /// </summary>
        /// <returns></returns>
        public string GetAaChange()
        {
            if (ReferenceAA == "" && AlternateAA == "")
            {
                if (CodonNum >= 0) { return "" + (CodonNum + 1); }
                return "";
            }

            if (ReferenceAA.Equals(AlternateAA)) { return AlternateAA + (CodonNum + 1); }
            return ReferenceAA + (CodonNum + 1) + AlternateAA;
        }

        /// <summary>
        /// Amino acid change string (old style)
        /// </summary>
        /// <returns></returns>
        public string GetAaChangeOld()
        {
            if (ReferenceAA == "" && AlternateAA == "") { return ""; }
            if (ReferenceAA.Equals(AlternateAA)) { return AlternateAA; }
            return (ReferenceAA == "" ? "-" : ReferenceAA) + "/" + (AlternateAA == "" ? "-" : AlternateAA);
        }

        /// <summary>
        /// CDS length (negative if there is none)
        /// </summary>
        /// <returns></returns>
        public long GetCdsLength()
        {
            Transcript tr = GetTranscript();
            return tr != null && tr.IsProteinCoding() ?
                tr.RetrieveCodingSequence().Count :
                -1;
        }

        /// <summary>
        /// Codon change string
        /// </summary>
        /// <returns></returns>
        public string GetCodonChange()
        {
            return CodonsRef == "" && CodonsAlt == "" ?
                "" :
                CodonsRef + "/" + CodonsAlt;
        }

        /// <summary>
        /// Codon change string (if it's not too long)
        /// </summary>
        /// <returns></returns>
        public string GetCodonChangeMax()
        {
            if (Variant.Length() > MAX_CODON_SEQUENCE_LEN) { return ""; }// Cap length in order not to make VCF files grow too much
            if (CodonsRef == "" && CodonsAlt == "") { return ""; }
            return CodonsRef + "/" + CodonsAlt;
        }

        /// <summary>
        /// Return impact of this effect
        /// </summary>
        /// <returns></returns>
        public EffectImpact GetEffectImpact()
        {
            //if ((variant != null) && (!variant.VarType == Variant.VariantType.()))
            //{
            //    // Not a change? => Modifier
            //    return EffectImpact.MODIFIER;
            //}
            return EffectImpacts.Min(); // Get effect's type highest impact
        }

        /// <summary>
        /// Highest effect type
        /// </summary>
        /// <returns></returns>
        public EffectType GetEffectType()
        {
            // Pick highest effect type
            return EffectTypes == null ? EffectType.NONE : EffectTypes.Min();
        }

        public string GetEffectTypeString(bool useFirstEffect)
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

        public string GetError()
        {
            return String.Join(",", Error);
        }

        /// <summary>
        /// Get exon (if any)
        /// </summary>
        /// <returns></returns>
        public Exon GetExon()
        {
            if (Marker != null)
            {
                if (Marker is Exon)
                {
                    return (Exon)Marker;
                }
                return (Exon)Marker.FindParent(typeof(Exon));
            }
            return null;
        }

        /// <summary>
        /// Return functional class of this effect (i.e.  NONSENSE, MISSENSE, SILENT or NONE)
        /// </summary>
        /// <returns></returns>
        public FunctionalClass GetFunctionalClass()
        {
            if (Variant.isSnv())
            {
                if (!AlternateAA.Equals(ReferenceAA))
                {
                    Translation.TranslateDnaCodon(CodonsAlt, out byte aa);
                    if (aa == Alphabets.Protein.Ter) { return FunctionalClass.NONSENSE; }

                    return FunctionalClass.MISSENSE;
                }
                if (!CodonsAlt.Equals(CodonsRef))
                {
                    return FunctionalClass.SILENT;
                }
            }

            return FunctionalClass.NONE;
        }

        /// <summary>
        /// Checks if this effect is nonsynonymous
        /// </summary>
        /// <returns></returns>
        public bool IsNonsynonymous()
        {
            return GetFunctionalClass().CompareTo(FunctionalClass.MISSENSE) >= 0;
        }

        public Gene GetGene()
        {
            if (Marker != null)
            {
                if (Marker is Gene) { return (Gene)Marker; }
                return (Gene)Marker.FindParent(typeof(Gene));
            }
            return null;
        }

        public string GetGeneRegion()
        {
            EffectType eff = GetGeneRegion(GetEffectType());
            if (eff == EffectType.TRANSCRIPT && IsExon()) { eff = EffectType.EXON; }
            return eff.ToString();
        }

        /// <summary>
        /// Get intron (if any)
        /// </summary>
        /// <returns></returns>
        public Intron GetIntron()
        {
            if (Marker != null)
            {
                if (Marker is Intron) { return (Intron)Marker; }
                return (Intron)Marker.FindParent(typeof(Intron));
            }
            return null;
        }

        public Interval GetMarker()
        {
            return Marker;
        }

        public Transcript GetTranscript()
        {
            if (Marker != null)
            {
                if (Marker is Transcript) { return (Transcript)Marker; }
                return (Transcript)Marker.FindParent(typeof(Transcript));
            }
            return null;
        }

        public Variant GetVariant()
        {
            return Variant;
        }

        public string GetWarning()
        {
            return String.Join(",", Warning);
        }

        public bool HasEffectImpact(EffectImpact effectImpact)
        {
            foreach (EffectImpact effImp in EffectImpacts)
            {
                if (effImp == effectImpact) { return true; }
            }
            return false;
        }

        public bool HasEffectType(EffectType effectType)
        {
            foreach (EffectType effType in EffectTypes)
            {
                if (effType == effectType) { return true; }
            }
            return false;
        }

        public bool IsExon()
        {
            return Marker is Exon || HasEffectType(EffectType.EXON_DELETED);
        }

        public bool IsIntergenic()
        {
            return HasEffectType(EffectType.INTERGENIC) || HasEffectType(EffectType.INTERGENIC_CONSERVED);
        }

        public bool IsIntron()
        {
            return HasEffectType(EffectType.INTRON) || HasEffectType(EffectType.INTRON_CONSERVED);
        }

        public bool IsMotif()
        {
            return HasEffectType(EffectType.MOTIF);
        }

        public bool IsMultipleGenes()
        {
            return false;
        }

        public bool IsNextProt()
        {
            return HasEffectType(EffectType.NEXT_PROT);
        }

        public bool IsRegulation()
        {
            return HasEffectType(EffectType.REGULATION);
        }

        public bool IsSpliceSite()
        {
            return HasEffectType(EffectType.SPLICE_SITE_DONOR) //
                    || HasEffectType(EffectType.SPLICE_SITE_ACCEPTOR) //
                    || HasEffectType(EffectType.SPLICE_SITE_REGION) //
                    || HasEffectType(EffectType.SPLICE_SITE_BRANCH) //
                    || HasEffectType(EffectType.SPLICE_SITE_BRANCH_U12); //
        }

        public bool IsSpliceSiteCore()
        {
            return HasEffectType(EffectType.SPLICE_SITE_DONOR) //
                    || HasEffectType(EffectType.SPLICE_SITE_ACCEPTOR); //
        }

        public bool IsSpliceSiteRegion()
        {
            return HasEffectType(EffectType.SPLICE_SITE_REGION);
        }

        public bool IsUtr3()
        {
            return HasEffectType(EffectType.UTR_3_PRIME) || HasEffectType(EffectType.UTR_3_DELETED);
        }

        public bool IsUtr5()
        {
            return HasEffectType(EffectType.UTR_5_PRIME) || HasEffectType(EffectType.UTR_5_DELETED);
        }

        public void Set(Interval marker, EffectType effectType, EffectImpact effectImpact, string message)
        {
            SetMarker(marker); // Use setter because it takes care of warnings
            SetEffectType(effectType);
            SetEffectImpact(effectImpact);
            this.Message = message;
        }

        /// <summary>
        /// Set codon change. Calculate effect type based on codon changes (for SNPs & MNPs)
        /// </summary>
        /// <param name="codonsOld"></param>
        /// <param name="codonsNew"></param>
        /// <param name="codonNum"></param>
        /// <param name="codonIndex"></param>
        public void SetCodons(string codonsOld, string codonsNew, int codonNum, int codonIndex)
        {
            CodonsRef = codonsOld;
            CodonsAlt = codonsNew;
            CodonNum = codonNum;
            CodonIndex = codonIndex;

            // Calculate amino acids
            if (codonsOld == "")
            {
                ReferenceAA = "";
            }
            else if (Translation.TranslateDnaCodon(codonsOld, out byte aa))
            {
                ReferenceAA = new string(new[] { (char)aa });
                //codonDegeneracy = codonTable.degenerate(codonsOld, codonIndex); // Calculate codon degeneracy
            }

            if (codonsNew == "")
            {
                AlternateAA = "";
            }
            else if (Translation.TranslateDnaCodon(codonsNew, out byte aa))
            {
                AlternateAA = new string(new[] { (char)aa });
            }
        }

        /// <summary>
        /// Set values for codons around change.
        /// </summary>
        /// <param name="codonsLeft"></param>
        /// <param name="codonsRight"></param>
        public void SetCodonsAround(string codonsLeft, string codonsRight)
        {
            CodonsAroundOld = codonsLeft.ToLower(CultureInfo.InvariantCulture) + CodonsRef.ToUpper(CultureInfo.InvariantCulture) + codonsRight.ToLower(CultureInfo.InvariantCulture);
            CodonsAroundNew = codonsLeft.ToLower(CultureInfo.InvariantCulture) + CodonsAlt.ToUpper(CultureInfo.InvariantCulture) + codonsRight.ToLower(CultureInfo.InvariantCulture);

            // Amino acids surrounding the ones changed
            string aasLeft = Translation.TranslateDnaCodon(codonsLeft, out byte aa1) ?
                new string(new char[] { (char)aa1 }) :
                "";
            string aasRigt = Translation.TranslateDnaCodon(codonsRight, out byte aa2) ?
                new string(new char[] { (char)aa2 }) :
                "";
            AroundOldAAs = aasLeft.ToLower(CultureInfo.InvariantCulture) + ReferenceAA.ToUpper(CultureInfo.InvariantCulture) + aasRigt.ToLower(CultureInfo.InvariantCulture);
            AroundNewAAs = aasLeft.ToLower(CultureInfo.InvariantCulture) + AlternateAA.ToUpper(CultureInfo.InvariantCulture) + aasRigt.ToLower(CultureInfo.InvariantCulture);
        }

        public void SetDistance(long distance)
        {
            this.Distance = distance;
        }

        /// <summary>
        /// Set effect using default impact
        /// </summary>
        /// <param name="effectType"></param>
        public void SetEffect(EffectType effectType)
        {
            SetEffectType(effectType);
            SetEffectImpact(EffectDictionary[effectType]);
        }

        public void SetEffectImpact(EffectImpact effectImpact)
        {
            EffectImpacts.Clear();
            EffectImpacts.Add(effectImpact);
        }

        public void SetEffectType(EffectType effectType)
        {
            EffectTypes.Clear();
            EffectTypes.Add(effectType);
        }

        /// <summary>
        /// Set marker. Add some warnings if the marker relates to incomplete transcripts
        /// </summary>
        /// <param name="marker"></param>
        public void SetMarker(Interval marker)
        {
            Marker = marker;

            Transcript transcript = GetTranscript();
            if (transcript != null)
            {
                // Transcript level errors or warnings
                AddErrorWarningInfo(transcript.sanityCheck(Variant));

                // Exon level errors or warnings
                Exon exon = GetExon();
                if (exon != null)
                {
                    AddErrorWarningInfo(exon.SanityCheck(Variant));
                }
            }
        }

        /// <summary>
        /// Add an error or warning
        /// </summary>
        /// <param name="errwarn"></param>
        public void AddErrorWarningInfo(ErrorWarningType errwarn)
        {
            if (errwarn == ErrorWarningType.NONE) return;

            if (errwarn.ToString().StartsWith("ERROR"))
            {
                Error.Add(errwarn);
            }
            else
            {
                Warning.Add(errwarn);
            }
        }

        public string ToStr()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Variant             : " + Variant.ToString());
            sb.Append("\n\tEffectTypes     : " + EffectTypes);
            sb.Append("\n\tEffectImpacts   : " + EffectImpacts);
            if (Marker != null) sb.Append("\n\tMarker          : " + Marker.ToString());
            if (CodonsRef != "" || CodonsAlt != "") sb.Append("\n\tCodons Ref/Alt  : " + CodonsRef + " / " + CodonsAlt);
            if (ReferenceAA != "" || AlternateAA != "") sb.Append("\n\tAA Ref/Alt      : " + ReferenceAA + " / " + AlternateAA);
            if (CDnaPos >= 0) sb.Append("\n\tcDnaPos         : " + CDnaPos);
            if (CodonNum >= 0) sb.Append("\n\tcodon num/index : " + CodonNum + " / " + CodonIndex);
            if (Error.Count > 0) sb.Append("\n\tError           : " + String.Join(",", Error));
            if (Warning.Count > 0) sb.Append("\n\tWarning         : " + String.Join(",", Warning));
            if (Message != "") sb.Append("\n\tMessage         : " + Message);

            return sb.ToString();
        }

        public override string ToString()
        {
            return ToString(false, false);
        }

        public string ToString(bool useSeqOntology, bool useHgvs)
        {
            // Get data to show
            string geneId = "", geneName = "", transcriptId = "", exonId = "", customId = "";
            string bioType = null;
            int exonRank = -1;

            if (Marker != null)
            {
                // Gene Id, name and biotype
                Gene gene = GetGene();
                Transcript tr = GetTranscript();

                // CDS size info
                if (gene != null)
                {
                    geneId = gene.ID;
                    //geneName = gene.getGeneName();
                    //bioType = (getBiotype() == null ? "" : getBiotype().ToString());
                }

                // Update trId
                if (tr != null) transcriptId = tr.ID;

                // Exon rank information
                Exon exon = GetExon();
                //if (exon != null)
                //{
                //    exonId = exon.ID;
                //    exonRank = exon.getRank();
                //}

                // Regulation
                //if (isRegulation()) bioType = ((Regulation)marker).getRegulationType();
            }

            // Add seqChage's ID
            //if (!variant.getId() == "") customId += variant.getId();

            // Add custom markers
            //if (marker != null && marker is Custom) customId += (customId == "" ? "" : ";") + marker.getId();

            // CDS length
            long cdsSize = GetCdsLength();

            string error = String.Join(",", Error);
            string warning = String.Join(",", Warning);
            string errWarn = error + (error == "" ? "" : "|") + warning;

            //string aaChange = "";

            //if (useHgvs) aaChange = getHgvs();
            //else aaChange = ((aaRef.Length + aaAlt.Length) > 0 ? aaRef + "/" + aaAlt : "");

            string aaChange = ((ReferenceAA.Length + AlternateAA.Length) > 0 ? ReferenceAA + "/" + AlternateAA : "");

            return errWarn //
                    + "\t" + geneId //
                    + "\t" + geneName //
                    + "\t" + bioType //
                    + "\t" + transcriptId //
                    + "\t" + exonId //
                    + "\t" + (exonRank >= 0 ? exonRank.ToString() : "") //
                    + "\t" + Effect(false, false, false, useSeqOntology, false) //
                    + "\t" + aaChange //
                    + "\t" + ((CodonsRef.Length + CodonsAlt.Length) > 0 ? CodonsRef + "/" + CodonsAlt : "") //
                    + "\t" + (CodonNum >= 0 ? (CodonNum + 1).ToString() : "") //
                                                                              //+ "\t" + (codonDegeneracy >= 0 ? codonDegeneracy + "" : "") //
                    + "\t" + (cdsSize >= 0 ? cdsSize.ToString() : "") //
                    + "\t" + (CodonsAroundOld.Length > 0 ? CodonsAroundOld + " / " + CodonsAroundNew : "") //
                    + "\t" + (AroundOldAAs.Length > 0 ? AroundOldAAs + " / " + AroundNewAAs : "") //
                    + "\t" + customId //
            ;
        }

        public string TranscriptAnnotation()
        {
            return "effect:" + GetFunctionalClass().ToString() + " " +
                    GetEffectType().ToString() + " " +
                    ReferenceAA + (CodonNum + 1).ToString() + AlternateAA + " " +
                    CodonsRef.ToLower(CultureInfo.InvariantCulture) + CodonNum.ToString() + CodonsAlt;
        }

        /// <summary>
        /// Get the simplest string describing the effect (this is mostly used for testcases)
        /// </summary>
        /// <param name="shortFormat"></param>
        /// <returns></returns>
        public string ToStringSimple(bool shortFormat)
        {
            string transcriptId = "";
            Transcript tr = GetTranscript();
            if (tr != null) transcriptId = tr.ID;

            string exonId = "";
            Exon exon = GetExon();
            //if (exon != null) exonId = exon.getId();

            string eff = Effect(shortFormat, true, true, false, false);
            if (eff != "") return eff;
            if (exonId != "") return exonId;
            if (transcriptId != "") return transcriptId;

            return "NO EFFECT";
        }

        private static int CompareNull(IComparable oThis, IComparable oThat)
        {
            return oThis != null && oThat != null ? oThis.CompareTo(oThat) :
                oThis == null && oThat == null ? 0 :
                oThis == null ? -1 : 1;
        }
    }
}