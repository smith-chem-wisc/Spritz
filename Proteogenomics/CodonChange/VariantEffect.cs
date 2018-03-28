using Bio;
using Bio.Algorithms.Translation;
using Bio.VCF;
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

        public Variant Variant { get; private set; }
        public List<EffectType> effectTypes { get; private set; }
        protected EffectType effectType { get; set; }
        protected List<EffectImpact> EffectImpacts { get; set; }
        protected EffectImpact effectImpact { get; set; }
        protected Interval marker { get; set; }
        protected HashSet<ErrorWarningType> error { get; set; } = new HashSet<ErrorWarningType>();
        protected HashSet<ErrorWarningType> warning { get; set; } = new HashSet<ErrorWarningType>();
        protected string message { get; set; } = ""; // Any message, warning or error?
        public string codonsRef { get; private set; } = ""; // Codon change information
        public string codonsAlt { get; private set; } = ""; // Codon change information
        protected string codonsAroundOld = "", codonsAroundNew = ""; // Codons around
        protected long distance = -1; // Distance metric
        public long cDnaPos { get; private set; } = -1; // Position in cDNA
        public int codonNum { get; private set; } = -1; // Codon number (negative number mens 'information not available')
        public int codonIndex { get; private set; } = -1; // Index within a codon (negative number mens 'information not available')
        public int codonDegeneracy { get; private set; } = -1; // Codon degeneracy (negative number mens 'information not available')
        public string aaRef { get; set; } = "";
        public string aaAlt { get; private set; } = ""; // Amino acid changes
        protected string aasAroundOld { get; set; } = "";
        protected string aasAroundNew { get; set; } = ""; // Amino acids around

        public VariantEffect(Variant variant)
        {
            Variant = variant;
            effectTypes = new List<EffectType>();
            EffectImpacts = new List<EffectImpact>();
        }

        public VariantEffect(Variant variant, Interval marker, EffectType effectType, EffectImpact effectImpact, string codonsOld, string codonsNew, int codonNum, int codonIndex, long cDnaPos)
        {
            Variant = variant;
            effectTypes = new List<EffectType>();
            EffectImpacts = new List<EffectImpact>();
            set(marker, effectType, effectImpact, "");
            setCodons(codonsOld, codonsNew, codonNum, codonIndex);
            this.cDnaPos = cDnaPos;
        }

        public void addEffect(EffectType effectType)
        {
            addEffectType(effectType);
            addEffectImpact(EffectDictionary[effectType]);
        }

        public void addEffectImpact(EffectImpact effectImpact)
        {
            EffectImpacts.Add(effectImpact);
            this.effectImpact = EffectImpact.MODIFIER;
        }

        public void addEffectType(EffectType effectType)
        {
            effectTypes.Add(effectType);
            this.effectType = EffectType.NONE;
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
        private string codonEffect(bool showAaChange, bool showBioType, bool useFirstEffect)
        {
            if ((marker == null) || (codonNum < 0)) { return ""; }

            if (!showAaChange) { return getEffectTypeString(useFirstEffect); }

            StringBuilder sb = new StringBuilder();
            sb.Append(getEffectTypeString(useFirstEffect));
            sb.Append("(");
            sb.Append(getAaChange());
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
            int comp = getEffectImpact().CompareTo(varEffOther.getEffectImpact());
            if (comp != 0) return comp;

            // Sort by effect
            comp = getEffectType().CompareTo(varEffOther.getEffectType());
            if (comp != 0) return comp;

            //---
            // Transcript based comparisons
            //---
            Transcript trThis = getTranscript();
            Transcript trOther = varEffOther.getTranscript();

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
            comp = compareNull(chrThis, chrOther);
            if (comp != 0) { return comp; }

            comp = startThis.CompareTo(startOther);
            if (comp != 0) { return comp; }

            comp = endThis.CompareTo(endOther);
            if (comp != 0) { return comp; }

            // Compare IDs
            comp = compareNull(idThis, idOther);
            if (comp != 0) { return comp; }

            //---
            // Marker based comparisons
            //---
            Interval mThis = getMarker();
            Interval mOther = varEffOther.getMarker();

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
            comp = compareNull(chrThis, chrOther);
            if (comp != 0) { return comp; }

            comp = startThis.CompareTo(startOther);
            if (comp != 0) { return comp; }

            comp = endThis.CompareTo(endOther);
            if (comp != 0) { return comp; }

            // Compare IDs
            comp = compareNull(idThis, idOther);
            if (comp != 0) { return comp; }

            //---
            // Variant based comparison
            //---
            return Variant.CompareTo(varEffOther.Variant);
        }

        private static int compareNull(IComparable oThis, IComparable oThat)
        {
            return oThis != null && oThat != null ? oThis.CompareTo(oThat) :
                oThis == null && oThat == null ? 0 :
                oThis == null ? -1 : 1;
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
        public string effect(bool shortFormat, bool showAaChange, bool showBioType, bool useSeqOntology, bool useFirstEffect)
        {
            string e = "";
            string codonEffect = this.codonEffect(showAaChange, showBioType, useFirstEffect); // Codon effect

            // Create effect string
            if (codonEffect != "") { e = codonEffect; }
            //else if (isRegulation()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + ((Regulation)marker).getName() + "]";
            //else if (isNextProt()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + VcfEffect.vcfEffSafe(((NextProt)marker).getId()) + "]"; // Make sure this 'id' is not dangerous in a VCF 'EFF' field
            //else if (isMotif()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + ((Motif)marker).getPwmId() + ":" + ((Motif)marker).getPwmName() + "]";
            //else if (isCustom())
            //{
            //    // Custom interval
            //    string label = ((Custom)marker).getLabel();
            //    double score = ((Custom)marker).getScore();
            //    if (score != double.NaN) label = label + ":" + score;
            //    if (label != "") label = "[" + label + "]";
            //    return getEffectTypestring(useSeqOntology, useFirstEffect) + label;
            //}
            else if (isIntergenic() || isIntron() || isSpliceSite()) { e = getEffectTypeString(useFirstEffect); }
            else if (message != "") { e = getEffectTypeString(useFirstEffect) + ": " + message; }
            else if (marker == null) { e = getEffectTypeString(useFirstEffect); }// There are cases when no marker is associated (e.g. "Out of chromosome", "No such chromosome", etc.)
            else { e = getEffectTypeString(useFirstEffect); }// + ": " + marker.getId();

            if (shortFormat) { e = e.Split(':')[0]; }

            return e;
        }

        /// <summary>
        /// Amino acid change string (HGVS style)
        /// </summary>
        /// <returns></returns>
        public string getAaChange()
        {
            if (aaRef == "" && aaAlt == "")
            {
                if (codonNum >= 0) { return "" + (codonNum + 1); }
                return "";
            }

            if (aaRef.Equals(aaAlt)) { return aaAlt + (codonNum + 1); }
            return aaRef + (codonNum + 1) + aaAlt;
        }

        /// <summary>
        /// Amino acid change string (old style)
        /// </summary>
        /// <returns></returns>
        public string getAaChangeOld()
        {
            if (aaRef == "" && aaAlt == "") { return ""; }
            if (aaRef.Equals(aaAlt)) { return aaAlt; }
            return (aaRef == "" ? "-" : aaRef) + "/" + (aaAlt == "" ? "-" : aaAlt);
        }

        /// <summary>
        /// Amino acid length (negative if there is none)
        /// @return Amino acid length (CDS length / 3 ) or '-1' if there is no CDS length
        /// </summary>
        /// <returns></returns>
        public int getAaLength()
        {
            int cdsLen = (int)getCdsLength();
            if (cdsLen < 0) return -1;

            int lenNoStop = Math.Max(0, cdsLen - 3); // Do not include the STOP codon
            return lenNoStop / 3;
        }

        /// <summary>
        /// Net AA change (InDels)
        /// </summary>
        /// <returns></returns>
        public string getAaNetChange()
        {
            string aaShort = aaRef.ToUpper(CultureInfo.InvariantCulture);
            string aaLong = aaAlt.ToUpper(CultureInfo.InvariantCulture);

            if (aaLong.Length < aaShort.Length)
            {
                string tmp = aaShort;
                aaShort = aaLong;
                aaLong = tmp;
            }

            if (aaLong.StartsWith(aaShort)) { return aaLong.Substring(aaShort.Length); }
            if (aaLong.EndsWith(aaLong)) { return aaLong.Substring(0, aaLong.Length - aaShort.Length); }
            if (aaShort == "") { return aaLong; }

            // Assumptions broken (may be this is not an InDel).
            return null;
        }

        /// <summary>
        /// Get biotype
        /// </summary>
        /// <returns></returns>
        //public BioType getBiotype()
        //{
        //    Gene gene = getGene();
        //    if (gene == null) return null;

        //    Transcript tr = getTranscript();
        //    if (tr != null) return tr.getBioType();
        //    else if (gene.getGenome().hasCodingInfo()) return BioType.coding(gene.isProteinCoding());

        //    return null;
        //}

        /// <summary>
        /// CDS length (negative if there is none)
        /// </summary>
        /// <returns></returns>
        public long getCdsLength()
        {
            // CDS size info
            Transcript tr = getTranscript();
            if ((tr != null) && tr.isProteinCoding()) { return tr.RetrieveCodingSequence().Count; }
            return -1;
        }

        /// <summary>
        /// Codon change string
        /// </summary>
        /// <returns></returns>
        public string getCodonChange()
        {
            if (codonsRef == "" && codonsAlt == "") { return ""; }
            return codonsRef + "/" + codonsAlt;
        }

        /// <summary>
        /// Codon change string (if it's not too long)
        /// </summary>
        /// <returns></returns>
        public string getCodonChangeMax()
        {
            if (Variant.Length() > MAX_CODON_SEQUENCE_LEN) { return ""; }// Cap length in order not to make VCF files grow too much
            if (codonsRef == "" && codonsAlt == "") { return ""; }
            return codonsRef + "/" + codonsAlt;
        }

        /// <summary>
        /// Return impact of this effect
        /// </summary>
        /// <returns></returns>
        public EffectImpact getEffectImpact()
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
        public EffectType getEffectType()
        {
            if (effectType != EffectType.NONE)
            {
                return effectType;
            }
            if (effectTypes == null)
            {
                return EffectType.NONE;
            }

            // Pick highest effect type
            effectType = EffectType.NONE;
            return effectTypes.Min();
        }

        public string getEffectTypeString(bool useFirstEffect)
        {
            if (effectTypes == null) { return ""; }

            // Show all effects
            StringBuilder sb = new StringBuilder();
            effectTypes.Sort();

            // More than one effect? Check for repeats
            HashSet<string> added = ((effectTypes.Count > 1) && (!useFirstEffect) ? new HashSet<string>() : null);

            // Create string
            foreach (EffectType et in effectTypes)
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

        public string getError()
        {
            return String.Join(",", error);
        }

        /// <summary>
        /// Get exon (if any)
        /// </summary>
        /// <returns></returns>
        public Exon getExon()
        {
            if (marker != null)
            {
                if (marker is Exon)
                {
                    return (Exon)marker;
                }
                return (Exon)marker.FindParent(typeof(Exon));
            }
            return null;
        }

        /// <summary>
        /// Return functional class of this effect (i.e.  NONSENSE, MISSENSE, SILENT or NONE)
        /// </summary>
        /// <returns></returns>
        public FunctionalClass getFunctionalClass()
        {
            if (Variant.isSnv())
            {
                if (!aaAlt.Equals(aaRef))
                {
                    Codons.TryLookup((byte)codonsAlt[0], (byte)codonsAlt[1], (byte)codonsAlt[2], out byte aa);
                    if (aa == Alphabets.Protein.Ter) { return FunctionalClass.NONSENSE; }

                    return FunctionalClass.MISSENSE;
                }
                if (!codonsAlt.Equals(codonsRef)) { return FunctionalClass.SILENT; }
            }

            return FunctionalClass.NONE;
        }

        public Gene getGene()
        {
            if (marker != null)
            {
                if (marker is Gene) { return (Gene)marker; }
                return (Gene)marker.FindParent(typeof(Gene));
            }
            return null;
        }

        public string getGeneRegion()
        {
            EffectType eff = getGeneRegion(getEffectType());
            if (eff == EffectType.TRANSCRIPT && isExon()) { eff = EffectType.EXON; }
            return eff.ToString();
        }

        public static EffectType getGeneRegion(EffectType type)
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

        public static List<Gene> getGenes()
        {
            return null;
        }

        /// <summary>
        /// Get genotype string
        /// </summary>
        /// <returns></returns>
        public Genotype getGenotype()
        {
            return Variant.Genotype;
        }

        ///**
        // * Change in HGVS notation
        // */
        //public string getHgvs()
        //{
        //    if (!Config.get().isHgvs()) return "";

        //    // Calculate protein level and dna level changes
        //    string hgvsProt = getHgvsProt();
        //    string hgvsDna = getHgvsDna();

        //    // Build output
        //    StringBuilder hgvs = new StringBuilder();
        //    if (hgvsProt != null) hgvs.Append(hgvsProt);
        //    if (hgvsDna != null)
        //    {
        //        if (hgvs.Length > 0) hgvs.Append('/');
        //        hgvs.Append(hgvsDna);
        //    }

        //    return hgvs.ToString();
        //}

        ///**
        // * Change in HGVS (Dna) notation
        // */
        //public string getHgvsDna()
        //{
        //    if (!Config.get().isHgvs()) return "";

        //    HgvsDna hgvsDna = new HgvsDna(this);
        //    string hgvs = hgvsDna.ToString();
        //    return hgvs != null ? hgvs : "";
        //}

        ///**
        // * Change in HGVS (Protein) notation
        // */
        //public string getHgvsProt()
        //{
        //    if (!Config.get().isHgvs()) return "";

        //    HgvsProtein hgvsProtein = new HgvsProtein(this);
        //    string hgvs = hgvsProtein.ToString();
        //    return hgvs != null ? hgvs : "";
        //}

        /**
         * Get intron (if any)
         */

        public Intron getIntron()
        {
            if (marker != null)
            {
                if (marker is Intron) { return (Intron)marker; }
                return (Intron)marker.FindParent(typeof(Intron));
            }
            return null;
        }

        public Interval getMarker()
        {
            return marker;
        }

        public Transcript getTranscript()
        {
            if (marker != null)
            {
                if (marker is Transcript) { return (Transcript)marker; }
                return (Transcript)marker.FindParent(typeof(Transcript));
            }
            return null;
        }

        public Variant getVariant()
        {
            return Variant;
        }

        public string getWarning()
        {
            return String.Join(",", warning);
        }

        /// <summary>
        /// Do we have an associated marker with additional annotations?
        /// </summary>
        /// <returns></returns>
        public bool hasAdditionalAnnotations()
        {
            //return getMarker() != null // Do we have a marker?
            //        && (getMarker() is Custom) // Is it 'custom'?
            //            && ((Custom)getMarker()).hasAnnotations(); // Does it have additional annotations?
            return false;
        }

        public bool hasEffectImpact(EffectImpact effectImpact)
        {
            foreach (EffectImpact effImp in EffectImpacts)
            {
                if (effImp == effectImpact) { return true; }
            }
            return false;
        }

        public bool hasEffectType(EffectType effectType)
        {
            foreach (EffectType effType in effectTypes)
            {
                if (effType == effectType) { return true; }
            }
            return false;
        }

        public bool hasError()
        {
            return error.Count > 0;
        }

        public bool hasWarning()
        {
            return warning.Count > 0;
        }

        public bool isCustom()
        {
            return getEffectType() == EffectType.CUSTOM;
        }

        public bool isExon()
        {
            return marker is Exon || hasEffectType(EffectType.EXON_DELETED);
        }

        public bool isIntergenic()
        {
            return hasEffectType(EffectType.INTERGENIC) || hasEffectType(EffectType.INTERGENIC_CONSERVED);
        }

        public bool isIntron()
        {
            return hasEffectType(EffectType.INTRON) || hasEffectType(EffectType.INTRON_CONSERVED);
        }

        public bool isMotif()
        {
            return hasEffectType(EffectType.MOTIF);
        }

        public bool isMultipleGenes()
        {
            return false;
        }

        public bool isNextProt()
        {
            return hasEffectType(EffectType.NEXT_PROT);
        }

        public bool isRegulation()
        {
            return hasEffectType(EffectType.REGULATION);
        }

        public bool isSpliceSite()
        {
            return hasEffectType(EffectType.SPLICE_SITE_DONOR) //
                    || hasEffectType(EffectType.SPLICE_SITE_ACCEPTOR) //
                    || hasEffectType(EffectType.SPLICE_SITE_REGION) //
                    || hasEffectType(EffectType.SPLICE_SITE_BRANCH) //
                    || hasEffectType(EffectType.SPLICE_SITE_BRANCH_U12); //
        }

        public bool isSpliceSiteCore()
        {
            return hasEffectType(EffectType.SPLICE_SITE_DONOR) //
                    || hasEffectType(EffectType.SPLICE_SITE_ACCEPTOR); //
        }

        public bool isSpliceSiteRegion()
        {
            return hasEffectType(EffectType.SPLICE_SITE_REGION);
        }

        public bool isUtr3()
        {
            return hasEffectType(EffectType.UTR_3_PRIME) || hasEffectType(EffectType.UTR_3_DELETED);
        }

        public bool isUtr5()
        {
            return hasEffectType(EffectType.UTR_5_PRIME) || hasEffectType(EffectType.UTR_5_DELETED);
        }

        public void set(Interval marker, EffectType effectType, EffectImpact effectImpact, string message)
        {
            setMarker(marker); // Use setter because it takes care of warnings
            setEffectType(effectType);
            setEffectImpact(effectImpact);
            this.message = message;
        }

        /// <summary>
        /// Set codon change. Calculate effect type based on codon changes (for SNPs & MNPs)
        /// </summary>
        /// <param name="codonsOld"></param>
        /// <param name="codonsNew"></param>
        /// <param name="codonNum"></param>
        /// <param name="codonIndex"></param>
        public void setCodons(string codonsOld, string codonsNew, int codonNum, int codonIndex)
        {
            codonsRef = codonsOld;
            codonsAlt = codonsNew;
            this.codonNum = codonNum;
            this.codonIndex = codonIndex;

            // Calculate amino acids
            if (codonsOld == "")
            {
                aaRef = "";
            }
            else if (Codons.TryLookup((byte)codonsOld[0], (byte)codonsOld[0], (byte)codonsOld[0], out byte aa))
            {
                aaRef = new string(new[] { (char)aa });
                //codonDegeneracy = codonTable.degenerate(codonsOld, codonIndex); // Calculate codon degeneracy
            }

            if (codonsNew == "")
            {
                aaAlt = "";
            }
            else if (Codons.TryLookup((byte)codonsNew[0], (byte)codonsNew[0], (byte)codonsNew[0], out byte aa))
            {
                aaAlt = new string(new[] { (char)aa });
            }
        }

        /// <summary>
        /// Set values for codons around change.
        /// </summary>
        /// <param name="codonsLeft"></param>
        /// <param name="codonsRight"></param>
        public void setCodonsAround(string codonsLeft, string codonsRight)
        {
            codonsAroundOld = codonsLeft.ToLower(CultureInfo.InvariantCulture) + codonsRef.ToUpper(CultureInfo.InvariantCulture) + codonsRight.ToLower(CultureInfo.InvariantCulture);
            codonsAroundNew = codonsLeft.ToLower(CultureInfo.InvariantCulture) + codonsAlt.ToUpper(CultureInfo.InvariantCulture) + codonsRight.ToLower(CultureInfo.InvariantCulture);

            // Amino acids surrounding the ones changed
            string aasLeft = Codons.TryLookup((byte)codonsLeft[0], (byte)codonsLeft[0], (byte)codonsLeft[0], out byte aa1) ?
                new string(new char[] { (char)aa1 }) :
                "";
            string aasRigt = Codons.TryLookup((byte)codonsRight[0], (byte)codonsRight[0], (byte)codonsRight[0], out byte aa2) ?
                new string(new char[] { (char)aa2 }) :
                "";
            aasAroundOld = aasLeft.ToLower(CultureInfo.InvariantCulture) + aaRef.ToUpper(CultureInfo.InvariantCulture) + aasRigt.ToLower(CultureInfo.InvariantCulture);
            aasAroundNew = aasLeft.ToLower(CultureInfo.InvariantCulture) + aaAlt.ToUpper(CultureInfo.InvariantCulture) + aasRigt.ToLower(CultureInfo.InvariantCulture);
        }

        public void setDistance(long distance)
        {
            this.distance = distance;
        }

        /// <summary>
        /// Set effect using default impact
        /// </summary>
        /// <param name="effectType"></param>
        public void setEffect(EffectType effectType)
        {
            setEffectType(effectType);
            setEffectImpact(EffectDictionary[effectType]);
        }

        public void setEffectImpact(EffectImpact effectImpact)
        {
            EffectImpacts.Clear();
            EffectImpacts.Add(effectImpact);
            this.effectImpact = EffectImpact.MODIFIER;
        }

        public void setEffectType(EffectType effectType)
        {
            effectTypes.Clear();
            effectTypes.Add(effectType);
            this.effectType = EffectType.NONE;
        }

        /// <summary>
        /// Set marker. Add some warnings if the marker relates to incomplete transcripts
        /// </summary>
        /// <param name="marker"></param>
        public void setMarker(Interval marker)
        {
            this.marker = marker;

            Transcript transcript = getTranscript();
            if (transcript != null)
            {
                // Transcript level errors or warnings
                addErrorWarningInfo(transcript.sanityCheck(Variant));

                // Exon level errors or warnings
                Exon exon = getExon();
                if (exon != null) addErrorWarningInfo(exon.SanityCheck(Variant));
            }
        }

        /// <summary>
        /// Add an error or warning
        /// </summary>
        /// <param name="errwarn"></param>
        public void addErrorWarningInfo(ErrorWarningType errwarn)
        {
            if (errwarn == ErrorWarningType.NONE) return;

            if (errwarn.ToString().StartsWith("ERROR"))
            {
                error.Add(errwarn);
            }
            else
            {
                warning.Add(errwarn);
            }
        }

        public string toStr()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("Variant             : " + Variant.ToString());
            sb.Append("\n\tEffectTypes     : " + effectTypes);
            sb.Append("\n\tEffectImpacts   : " + EffectImpacts);
            if (marker != null) sb.Append("\n\tMarker          : " + marker.ToString());
            if (codonsRef != "" || codonsAlt != "") sb.Append("\n\tCodons Ref/Alt  : " + codonsRef + " / " + codonsAlt);
            if (aaRef != "" || aaAlt != "") sb.Append("\n\tAA Ref/Alt      : " + aaRef + " / " + aaAlt);
            if (cDnaPos >= 0) sb.Append("\n\tcDnaPos         : " + cDnaPos);
            if (codonNum >= 0) sb.Append("\n\tcodon num/index : " + codonNum + " / " + codonIndex);
            if (error.Count > 0) sb.Append("\n\tError           : " + String.Join(",", error));
            if (warning.Count > 0) sb.Append("\n\tWarning         : " + String.Join(",", warning));
            if (message != "") sb.Append("\n\tMessage         : " + message);

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

            if (marker != null)
            {
                // Gene Id, name and biotype
                Gene gene = getGene();
                Transcript tr = getTranscript();

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
                Exon exon = getExon();
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
            long cdsSize = getCdsLength();

            string errWarn = error + (getError() == "" ? "" : "|") + warning;

            //string aaChange = "";

            //if (useHgvs) aaChange = getHgvs();
            //else aaChange = ((aaRef.Length + aaAlt.Length) > 0 ? aaRef + "/" + aaAlt : "");

            string aaChange = ((aaRef.Length + aaAlt.Length) > 0 ? aaRef + "/" + aaAlt : "");

            return errWarn //
                    + "\t" + geneId //
                    + "\t" + geneName //
                    + "\t" + bioType //
                    + "\t" + transcriptId //
                    + "\t" + exonId //
                    + "\t" + (exonRank >= 0 ? exonRank.ToString() : "") //
                    + "\t" + effect(false, false, false, useSeqOntology, false) //
                    + "\t" + aaChange //
                    + "\t" + ((codonsRef.Length + codonsAlt.Length) > 0 ? codonsRef + "/" + codonsAlt : "") //
                    + "\t" + (codonNum >= 0 ? (codonNum + 1).ToString() : "") //
                                                                              //+ "\t" + (codonDegeneracy >= 0 ? codonDegeneracy + "" : "") //
                    + "\t" + (cdsSize >= 0 ? cdsSize.ToString() : "") //
                    + "\t" + (codonsAroundOld.Length > 0 ? codonsAroundOld + " / " + codonsAroundNew : "") //
                    + "\t" + (aasAroundOld.Length > 0 ? aasAroundOld + " / " + aasAroundNew : "") //
                    + "\t" + customId //
            ;
        }

        /// <summary>
        /// Get the simplest string describing the effect (this is mostly used for testcases)
        /// </summary>
        /// <param name="shortFormat"></param>
        /// <returns></returns>
        public string ToStringSimple(bool shortFormat)
        {
            string transcriptId = "";
            Transcript tr = getTranscript();
            if (tr != null) transcriptId = tr.ID;

            string exonId = "";
            Exon exon = getExon();
            //if (exon != null) exonId = exon.getId();

            string eff = effect(shortFormat, true, true, false, false);
            if (eff != "") return eff;
            if (exonId != "") return exonId;
            if (transcriptId != "") return transcriptId;

            return "NO EFFECT";
        }
    }
}