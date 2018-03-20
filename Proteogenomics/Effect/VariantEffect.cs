using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Proteogenomics
{
    public class VariantEffect
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

        protected Variant variant { get; set; }
        protected List<EffectType> effectTypes { get; set; }
        protected EffectType effectType { get; set; }
        protected List<EffectImpact> effectImpacts { get; set; }
        protected EffectImpact effectImpact { get; set; }
        protected Interval marker { get; set; }
        protected string error { get; set; } = "";
        protected string warning { get; set; } = "";
        protected string message { get; set; } = ""; // Any message, warning or error?
        protected string codonsRef = "", codonsAlt = ""; // Codon change information
        protected string codonsAroundOld = "", codonsAroundNew = ""; // Codons around
        protected int distance = -1; // Distance metric
        protected int cDnaPos = -1; // Position in cDNA
        protected int codonNum = -1; // Codon number (negative number mens 'information not available')
        protected int codonIndex = -1; // Index within a codon (negative number mens 'information not available')
        protected int codonDegeneracy = -1; // Codon degeneracy (negative number mens 'information not available')
        protected string aaRef { get; set; } = "";
        protected string aaAlt { get; set; } = ""; // Amino acid changes
        protected string aasAroundOld { get; set; } = "";
        protected string aasAroundNew { get; set; } = ""; // Amino acids around

        public VariantEffect(Variant variant)
        {
            this.variant = variant;
            effectTypes = new List<EffectType>();
            effectImpacts = new List<EffectImpact>();
        }

        public VariantEffect(Variant variant, Interval marker, EffectType effectType, EffectImpact effectImpact, string codonsOld, string codonsNew, int codonNum, int codonIndex, int cDnaPos)
        {
            this.variant = variant;
            effectTypes = new List<EffectType>();
            effectImpacts = new List<EffectImpact>();
            set(marker, effectType, effectImpact, "");
            setCodons(codonsOld, codonsNew, codonNum, codonIndex);
            this.cDnaPos = cDnaPos;
        }

        public void addEffect(EffectType effectType)
        {
            addEffectType(effectType);
            addEffectImpact(effectType.effectImpact());
        }

        public void addEffectImpact(EffectImpact effectImpact)
        {
            effectImpacts.Add(effectImpact);
            this.effectImpact = null;
        }

        public void addEffectType(EffectType effectType)
        {
            effectTypes.Add(effectType);
            this.effectType = null;
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
        string codonEffect(bool showAaChange, bool showBioType, bool useSeqOntology, bool useFirstEffect)
        {
            if ((marker == null) || (codonNum < 0)) return "";

            if (!showAaChange) return getEffectTypestring(useSeqOntology, useFirstEffect);

            StringBuilder sb = new StringBuilder();
            sb.Append(getEffectTypestring(useSeqOntology, useFirstEffect));
            sb.Append("(");
            sb.Append(getAaChange());
            sb.Append(")");

            return sb.Tostring();
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
            string codonEffect = codonEffect(showAaChange, showBioType, useSeqOntology, useFirstEffect); // Codon effect

            // Create effect string
            if (!codonEffect.isEmpty()) e = codonEffect;
            else if (isRegulation()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + ((Regulation)marker).getName() + "]";
            else if (isNextProt()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + VcfEffect.vcfEffSafe(((NextProt)marker).getId()) + "]"; // Make sure this 'id' is not dangerous in a VCF 'EFF' field
            else if (isMotif()) return getEffectTypestring(useSeqOntology, useFirstEffect) + "[" + ((Motif)marker).getPwmId() + ":" + ((Motif)marker).getPwmName() + "]";
            else if (isCustom())
            {
                // Custom interval
                string label = ((Custom)marker).getLabel();
                double score = ((Custom)marker).getScore();
                if (!Double.isNaN(score)) label = label + ":" + score;
                if (!label.isEmpty()) label = "[" + label + "]";
                return getEffectTypestring(useSeqOntology, useFirstEffect) + label;
            }
            else if (isIntergenic() || isIntron() || isSpliceSite()) e = getEffectTypestring(useSeqOntology, useFirstEffect);
            else if (!message.isEmpty()) e = getEffectTypestring(useSeqOntology, useFirstEffect) + ": " + message;
            else if (marker == null) e = getEffectTypestring(useSeqOntology, useFirstEffect); // There are cases when no marker is associated (e.g. "Out of chromosome", "No such chromosome", etc.)
            else e = getEffectTypestring(useSeqOntology, useFirstEffect) + ": " + marker.getId();

            if (shortFormat) e = e.split(":")[0];

            return e;
        }

        public string getAaAlt()
        {
            return aaAlt;
        }

        /**
         * Amino acid change string (HGVS style)
         */
        public string getAaChange()
        {
            if (aaRef.isEmpty() && aaAlt.isEmpty())
            {
                if (codonNum >= 0) return "" + (codonNum + 1);
                return "";
            }

            if (aaRef.equals(aaAlt)) return aaAlt + (codonNum + 1);
            return aaRef + (codonNum + 1) + aaAlt;
        }

        /**
         * Amino acid change string (old style)
         */
        public string getAaChangeOld()
        {
            if (aaRef.isEmpty() && aaAlt.isEmpty()) return "";
            if (aaRef.equals(aaAlt)) return aaAlt;
            return (aaRef.isEmpty() ? "-" : aaRef) + "/" + (aaAlt.isEmpty() ? "-" : aaAlt);
        }

        /// <summary>
        /// Amino acid length (negative if there is none)
        /// @return Amino acid length (CDS length / 3 ) or '-1' if there is no CDS length
        /// </summary>
        /// <returns></returns>
        public int getAaLength()
        {
            int cdsLen = getCdsLength();
            if (cdsLen < 0) return -1;

            int lenNoStop = Math.max(0, cdsLen - 3); // Do not include the STOP codon
            return lenNoStop / 3;
        }

        /// <summary>
        /// Net AA change (InDels)
        /// </summary>
        /// <returns></returns>
        public string getAaNetChange()
        {
            string aaShort = getAaRef().toUpperCase();
            string aaLong = getAaAlt().toUpperCase();

            if (aaLong.length() < aaShort.length())
            {
                string tmp = aaShort;
                aaShort = aaLong;
                aaLong = tmp;
            }

            if (aaLong.startsWith(aaShort)) return aaLong.substring(aaShort.length());
            if (aaLong.endsWith(aaLong)) return aaLong.substring(0, aaLong.length() - aaShort.length());
            if (aaShort.isEmpty()) return aaLong;

            // Assumptions broken (may be this is not an InDel).
            return null;
        }

        public string getAaRef()
        {
            return aaRef;
        }

        /**
         * Get biotype
         */
        public BioType getBiotype()
        {
            Gene gene = getGene();
            if (gene == null) return null;

            Transcript tr = getTranscript();
            if (tr != null) return tr.getBioType();
            else if (gene.getGenome().hasCodingInfo()) return BioType.coding(gene.isProteinCoding());

            return null;
        }

        public int getcDnaPos()
        {
            return cDnaPos;
        }

        /// <summary>
        /// CDS length (negative if there is none)
        /// </summary>
        /// <returns></returns>
        public int getCdsLength()
        {
            // CDS size info
            Transcript tr = getTranscript();
            if ((tr != null) && tr.isProteinCoding()) return tr.cds().length();
            return -1;
        }

        /// <summary>
        /// Codon change string
        /// </summary>
        /// <returns></returns>
        public string getCodonChange()
        {
            if (codonsRef.isEmpty() && codonsAlt.isEmpty()) return "";
            return codonsRef + "/" + codonsAlt;
        }

        /// <summary>
        /// Codon change string (if it's not too long)
        /// </summary>
        /// <returns></returns>
        public string getCodonChangeMax()
        {
            if (variant.Length > MAX_CODON_SEQUENCE_LEN) return ""; // Cap length in order not to make VCF files grow too much
            if (codonsRef.isEmpty() && codonsAlt.isEmpty()) return "";
            return codonsRef + "/" + codonsAlt;
        }

        public int getCodonIndex()
        {
            return codonIndex;
        }

        public int getCodonNum()
        {
            return codonNum;
        }

        public string getCodonsAlt()
        {
            return codonsAlt;
        }

        public string getCodonsRef()
        {
            return codonsRef;
        }

        public int getDistance()
        {
            return distance;
        }

        /**
         * Return impact of this effect
         */
        public EffectImpact getEffectImpact()
        {
            if (effectImpact == null)
            {
                if ((variant != null) && (!variant.isVariant()))
                {
                    // Not a change? => Modifier
                    effectImpact = EffectImpact.MODIFIER;
                }
                else
                {
                    // Get efefct's type highest impact
                    effectImpact = EffectImpact.MODIFIER;
                    for (EffectImpact eimp : effectImpacts)
                        if (eimp.compareTo(effectImpact) < 0) effectImpact = eimp;
                }
            }

            return effectImpact;
        }

        /**
         * Highest effect type
         */
        public EffectType getEffectType()
        {
            if (effectType != null) return effectType;
            if (effectTypes == null || effectTypes.isEmpty()) return EffectType.NONE;

            // Pick highest effect type
            effectType = EffectType.NONE;
            for (EffectType et : effectTypes)
                if (et.compareTo(effectType) < 0) effectType = et;

            return effectType;
        }

        /**
         * Highest effect type
         */
        public List<EffectType> getEffectTypes()
        {
            return effectTypes;
        }

        public string getEffectTypestring(bool useSeqOntology)
        {
            return getEffectTypestring(useSeqOntology, false, EffFormatVersion.FORMAT_EFF_4);
        }

        public string getEffectTypestring(bool useSeqOntology, bool useFirstEffect)
        {
            return getEffectTypestring(useSeqOntology, useFirstEffect, EffFormatVersion.FORMAT_EFF_4);
        }

        /**
         * Get Effect Type as a string
         */
        public string getEffectTypestring(bool useSeqOntology, bool useFirstEffect, EffFormatVersion formatVersion)
        {
            if (effectTypes == null) return "";

            // Show all effects
            StringBuilder sb = new StringBuilder();
            Collections.sort(effectTypes);

            // More than one effect? Check for repeats
            Set<string> added = ((effectTypes.size() > 1) && (!useFirstEffect) ? new HashSet<string>() : null);

            // Create string
            for (EffectType et : effectTypes)
            {
                string eff = (useSeqOntology ? et.toSequenceOntology(formatVersion, variant) : et.tostring());

                // Make sure we don't add the same effect twice
                if (added == null || added.add(eff))
                {
                    if (sb.length() > 0) sb.append(formatVersion.separator());
                    sb.append(eff);
                }

                // Only use first effect?
                if (useFirstEffect) return sb.tostring();
            }

            return sb.tostring();
        }

        public string getError()
        {
            return error;
        }

        /**
         * Get exon (if any)
         */
        public Exon getExon()
        {
            if (marker != null)
            {
                if (marker is Exon)
                {
                    return (Exon)marker;
                }
                return (Exon)marker.findParent(Exon);
            }
		    return null;
        }

        /**
         * Return functional class of this effect (i.e.  NONSENSE, MISSENSE, SILENT or NONE)
         */
        public FunctionalClass getFunctionalClass()
        {
            if (variant.isSnp())
            {
                if (!aaAlt.equals(aaRef))
                {
                    CodonTable codonTable = marker.codonTable();
                    if (codonTable.isStop(codonsAlt)) return FunctionalClass.NONSENSE;

                    return FunctionalClass.MISSENSE;
                }
                if (!codonsAlt.equals(codonsRef)) return FunctionalClass.SILENT;
            }

            return FunctionalClass.NONE;
        }

        public Gene getGene()
        {
            if (marker != null)
            {
                if (marker instanceof Gene) return (Gene)marker;
                return (Gene)marker.findParent(Gene.class);
		        }
		        return null;
	        }

	        public string getGeneRegion()
        {
            EffectType eff = getEffectType().getGeneRegion();
            if (eff == EffectType.TRANSCRIPT && isExon()) eff = EffectType.EXON;
            return eff.tostring();
        }

        public List<Gene> getGenes()
        {
            return null;
        }

        /**
         * Get genotype string
         */
        public string getGenotype()
        {
            if (variant == null) return "";
            return variant.getGenotype();
        }

        /**
         * Change in HGVS notation
         */
        public string getHgvs()
        {
            if (!Config.get().isHgvs()) return "";

            // Calculate protein level and dna level changes
            string hgvsProt = getHgvsProt();
            string hgvsDna = getHgvsDna();

            // Build output
            StringBuilder hgvs = new StringBuilder();
            if (hgvsProt != null) hgvs.append(hgvsProt);
            if (hgvsDna != null)
            {
                if (hgvs.length() > 0) hgvs.append('/');
                hgvs.append(hgvsDna);
            }

            return hgvs.tostring();
        }

        /**
         * Change in HGVS (Dna) notation
         */
        public string getHgvsDna()
        {
            if (!Config.get().isHgvs()) return "";

            HgvsDna hgvsDna = new HgvsDna(this);
            string hgvs = hgvsDna.tostring();
            return hgvs != null ? hgvs : "";
        }

        /**
         * Change in HGVS (Protein) notation
         */
        public string getHgvsProt()
        {
            if (!Config.get().isHgvs()) return "";

            HgvsProtein hgvsProtein = new HgvsProtein(this);
            string hgvs = hgvsProtein.tostring();
            return hgvs != null ? hgvs : "";
        }

        /**
         * Get intron (if any)
         */
        public Intron getIntron()
        {
            if (marker != null)
            {
                if (marker instanceof Intron) return (Intron)marker;
                return (Intron)marker.findParent(Intron.class);
		        }
		        return null;
	        }

	        public Marker getMarker()
        {
            return marker;
        }

        public Transcript getTranscript()
        {
            if (marker != null)
            {
                if (marker instanceof Transcript) return (Transcript)marker;
                return (Transcript)marker.findParent(Transcript.class);
		        }
		        return null;
	        }

	        public Variant getVariant()
        {
            return variant;
        }

        public string getWarning()
        {
            return warning;
        }

        /**
         * Do we have an associated marker with additional annotations?
         */
        public bool hasAdditionalAnnotations()
        {
            return getMarker() != null // Do we have a marker?
                    && (getMarker() instanceof Custom) // Is it 'custom'?

                        && ((Custom)getMarker()).hasAnnotations() // Does it have additional annotations?
                ;
        }

        public bool hasEffectImpact(EffectImpact effectImpact)
        {
            for (EffectImpact effImp : effectImpacts)
                if (effImp == effectImpact) return true;
            return false;
        }

        public bool hasEffectType(EffectType effectType)
        {
            for (EffectType effType : effectTypes)
                if (effType == effectType) return true;
            return false;
        }

        public bool hasError()
        {
            return (error != null) && (!error.isEmpty());
        }

        public bool hasWarning()
        {
            return (warning != null) && (!warning.isEmpty());
        }

        public bool isCustom()
        {
            return getEffectType() == EffectType.CUSTOM;
        }

        public bool isExon()
        {
            return (marker instanceof Exon) || hasEffectType(EffectType.EXON_DELETED);
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
                    || hasEffectType(EffectType.SPLICE_SITE_BRANCH_U12) //
            ;
        }

        public bool isSpliceSiteCore()
        {
            return hasEffectType(EffectType.SPLICE_SITE_DONOR) //
                    || hasEffectType(EffectType.SPLICE_SITE_ACCEPTOR) //
            ;
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

        public void set(Marker marker, EffectType effectType, EffectImpact effectImpact, string message)
        {
            setMarker(marker); // Use setter because it takes care of warnings
            setEffectType(effectType);
            setEffectImpact(effectImpact);
            this.message = message;
        }

        /**
         * Set codon change. Calculate effect type based on codon changes (for SNPs & MNPs)
         */
        public void setCodons(string codonsOld, string codonsNew, int codonNum, int codonIndex)
        {
            codonsRef = codonsOld;
            codonsAlt = codonsNew;
            this.codonNum = codonNum;
            this.codonIndex = codonIndex;

            CodonTable codonTable = marker.codonTable();

            // Calculate amino acids
            if (codonsOld.isEmpty()) aaRef = "";
            else
            {
                aaRef = codonTable.aa(codonsOld);
                codonDegeneracy = codonTable.degenerate(codonsOld, codonIndex); // Calculate codon degeneracy
            }

            if (codonsNew.isEmpty()) aaAlt = "";
            else aaAlt = codonTable.aa(codonsNew);
        }

        /**
         * Set values for codons around change.
         */
        public void setCodonsAround(string codonsLeft, string codonsRight)
        {
            codonsAroundOld = codonsLeft.toLowerCase() + codonsRef.toUpperCase() + codonsRight.toLowerCase();
            codonsAroundNew = codonsLeft.toLowerCase() + codonsAlt.toUpperCase() + codonsRight.toLowerCase();

            // Amino acids surrounding the ones changed
            CodonTable codonTable = marker.codonTable();
            string aasLeft = codonTable.aa(codonsLeft);
            string aasRigt = codonTable.aa(codonsRight);
            aasAroundOld = aasLeft.toLowerCase() + aaRef.toUpperCase() + aasRigt.toLowerCase();
            aasAroundNew = aasLeft.toLowerCase() + aaAlt.toUpperCase() + aasRigt.toLowerCase();
        }

        public void setDistance(int distance)
        {
            this.distance = distance;
        }

        /**
         * Set effect using default impact
         */
        public void setEffect(EffectType effectType)
        {
            setEffectType(effectType);
            setEffectImpact(effectType.effectImpact());
        }

        public void setEffectImpact(EffectImpact effectImpact)
        {
            effectImpacts.clear();
            effectImpacts.add(effectImpact);
            this.effectImpact = null;
        }

        public void setEffectType(EffectType effectType)
        {
            effectTypes.clear();
            effectTypes.add(effectType);
            this.effectType = null;
        }

        /**
         * Set marker. Add some warnings if the marker relates to incomplete transcripts
         */
        public void setMarker(Marker marker)
        {
            this.marker = marker;

            Transcript transcript = getTranscript();
            if (transcript != null)
            {
                // Transcript level errors or warnings
                addErrorWarningInfo(transcript.sanityCheck(variant));

                // Exon level errors or warnings
                Exon exon = getExon();
                if (exon != null) addErrorWarningInfo(exon.sanityCheck(variant));
            }
        }

        public string toStr()
        {
            StringBuilder sb = new StringBuilder();
            sb.append("Variant             : " + variant.toStr());
            sb.append("\n\tEffectTypes     : " + effectTypes);
            sb.append("\n\tEffectImpacts   : " + effectImpacts);
            if (marker != null) sb.append("\n\tMarker          : " + marker.toStr());
            if (!codonsRef.isEmpty() || !codonsAlt.isEmpty()) sb.append("\n\tCodons Ref/Alt  : " + codonsRef + " / " + codonsAlt);
            if (!aaRef.isEmpty() || !aaAlt.isEmpty()) sb.append("\n\tAA Ref/Alt      : " + aaRef + " / " + aaAlt);
            if (cDnaPos >= 0) sb.append("\n\tcDnaPos         : " + cDnaPos);
            if (codonNum >= 0) sb.append("\n\tcodon num/index : " + codonNum + " / " + codonIndex);
            if (!error.isEmpty()) sb.append("\n\tError           : " + error);
            if (!warning.isEmpty()) sb.append("\n\tWarning         : " + warning);
            if (!message.isEmpty()) sb.append("\n\tMessage         : " + message);

            return sb.tostring();
        }

        @Override
            public string tostring()
        {
            return tostring(false, false);
        }

        public string tostring(bool useSeqOntology, bool useHgvs)
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
                    geneId = gene.getId();
                    geneName = gene.getGeneName();
                    bioType = (getBiotype() == null ? "" : getBiotype().tostring());
                }

                // Update trId
                if (tr != null) transcriptId = tr.getId();

                // Exon rank information
                Exon exon = getExon();
                if (exon != null)
                {
                    exonId = exon.getId();
                    exonRank = exon.getRank();
                }

                // Regulation
                if (isRegulation()) bioType = ((Regulation)marker).getRegulationType();
            }

            // Add seqChage's ID
            if (!variant.getId().isEmpty()) customId += variant.getId();

            // Add custom markers
            if ((marker != null) && (marker instanceof Custom)) customId += (customId.isEmpty() ? "" : ";") + marker.getId();

            // CDS length
            int cdsSize = getCdsLength();

            string errWarn = error + (error.isEmpty() ? "" : "|") + warning;

            string aaChange = "";
            if (useHgvs) aaChange = getHgvs();
            else aaChange = ((aaRef.length() + aaAlt.length()) > 0 ? aaRef + "/" + aaAlt : "");

            return errWarn //
                    + "\t" + geneId //
                    + "\t" + geneName //
                    + "\t" + bioType //
                    + "\t" + transcriptId //
                    + "\t" + exonId //
                    + "\t" + (exonRank >= 0 ? exonRank : "") //
                    + "\t" + effect(false, false, false, useSeqOntology, false) //
                    + "\t" + aaChange //
                    + "\t" + ((codonsRef.length() + codonsAlt.length()) > 0 ? codonsRef + "/" + codonsAlt : "") //
                    + "\t" + (codonNum >= 0 ? (codonNum + 1) : "") //
                    + "\t" + (codonDegeneracy >= 0 ? codonDegeneracy + "" : "") //
                    + "\t" + (cdsSize >= 0 ? cdsSize : "") //
                    + "\t" + (codonsAroundOld.length() > 0 ? codonsAroundOld + " / " + codonsAroundNew : "") //
                    + "\t" + (aasAroundOld.length() > 0 ? aasAroundOld + " / " + aasAroundNew : "") //
                    + "\t" + customId //
            ;
        }

        /**
         * Get the simplest string describing the effect (this is mostly used for testcases)
         */
        public string tostringSimple(bool shortFormat)
        {
            string transcriptId = "";
            Transcript tr = getTranscript();
            if (tr != null) transcriptId = tr.getId();

            string exonId = "";
            Exon exon = getExon();
            if (exon != null) exonId = exon.getId();

            string eff = effect(shortFormat, true, true, false, false);
            if (!eff.isEmpty()) return eff;
            if (!exonId.isEmpty()) return exonId;
            if (!transcriptId.isEmpty()) return transcriptId;

            return "NO EFFECT";
        }
    }
}
