using Bio;
using Bio.Extensions;
using Bio.VCF;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class Variant
        : Interval
    {
        #region Public Enum

        public enum VariantType
        {
            SNV // Single nucleotide variant (i.e. 1 base is changed)
            , MNP // Multiple nucleotide polymorphism (i.e. several bases are changed)
            , INS // Insertion (i.e. some bases added)
            , DEL // Deletion (some bases removed)
            , MIXED // A mixture of insertion, deletions, SNPs and or MNPs (a.k.a. substitution)
            , INV // Inversion (structural variant)
            , DUP // Duplication (structural variant)
            , BND // Break-ends (rearrangement)
            , INTERVAL // Just analyze interval hits. Not a variant (e.g. BED input format)
        }

        #endregion Public Enum

        #region Public Properties

        /// <summary>
        /// Original variant context
        /// </summary>
        public VariantContext VariantContext { get; set; }

        /// <summary>
        /// Type of this variant
        /// </summary>
        public VariantType VarType { get; set; }

        /// <summary>
        /// The reference allele
        /// </summary>
        public Allele ReferenceAllele { get; set; }

        /// <summary>
        /// The alternate allele
        /// </summary>
        public Allele AlternateAllele { get; set; }

        /// <summary>
        /// The reference allele as a string
        /// </summary>
        public string ReferenceAlleleString { get; set; }

        /// <summary>
        /// The alternate allele as a string
        /// </summary>
        public string AlternateAlleleString { get; set; }

        /// <summary>
        /// The genotype of this variation
        /// </summary>
        public Genotype Genotype { get; set; }

        /// <summary>
        /// The type of variation: heterozygous/homozygous
        /// </summary>
        public GenotypeType GenotypeType { get; set; }

        /// <summary>
        /// Depth across this position
        /// </summary>
        public int ReadDepth { get; set; }

        /// <summary>
        /// Depth of reads across this position containing the reference allele
        /// </summary>
        public int ReferenceAlleleDepth { get; set; }

        /// <summary>
        /// Depth of reads across this position containing the alternate allele
        /// </summary>
        public int AlternateAlleleDepth { get; set; }

        /// <summary>
        /// Text annotation for protein entries
        /// </summary>
        public string Annotation { get; set; }

        /// <summary>
        /// SnpEff annotaitons for this location.
        /// </summary>
        public List<SnpEffAnnotation> SnpEffAnnotations { get; } = new List<SnpEffAnnotation>();

        #endregion Public Properties

        #region Public Constructor

        public Variant(VariantContext variant, int zeroBasedAlleleIndex)
            : base(variant.Chr, "+", variant.Start, variant.End)
        {
            ReferenceAllele = variant.Reference;
            ReferenceAlleleString = variant.Reference.BaseString;
            Genotype = variant.Genotypes.First(); // assumes just one sample
            GenotypeType = Genotype.Type;
            ReadDepth = Genotype.DP;
            ReferenceAlleleDepth = Genotype.AD[0];
            AlternateAlleleDepth = Genotype.AD[1];
            AlternateAllele = Genotype.GetAllele(1);
            AlternateAlleleString = variant.AlternateAlleles[zeroBasedAlleleIndex].BaseString;
            CalculateType();

            if (variant.Attributes.TryGetValue("ANN", out object snpEffAnnotationObj))
            {
                if (snpEffAnnotationObj as string[] != null) SnpEffAnnotations = (snpEffAnnotationObj as string[]).Select(x => new SnpEffAnnotation(this, x)).ToList();
                if (snpEffAnnotationObj as string != null) SnpEffAnnotations = new List<SnpEffAnnotation> { new SnpEffAnnotation(this, snpEffAnnotationObj as string) };
            }
        }

        #endregion Public Constructor

        #region Public Methods

        public void CalculateType()
        {
            // Sanity check
            if (AlternateAlleleString.IndexOf(',') >= 0 || AlternateAlleleString.IndexOf('/') >= 0)
                throw new Exception("Variants with multiple ALTs are not allowed (ALT: '" + AlternateAlleleString + "')");

            // Remove leading char (we still have some test cases using old TXT format)
            if (ReferenceAlleleString == "*") ReferenceAlleleString = "";

            // Possibly some old formatting with +, -, and =
            if (AlternateAlleleString.StartsWith("+"))
            {
                // Insertion
                AlternateAlleleString = ReferenceAlleleString + AlternateAlleleString.Substring(1);
            }
            else if (AlternateAlleleString.StartsWith("-"))
            {
                // Deletion
                ReferenceAlleleString = AlternateAlleleString.Substring(1);
                AlternateAlleleString = "";
            }
            else if (AlternateAlleleString.StartsWith("="))
            {
                // Mixed variant
                AlternateAlleleString = AlternateAlleleString.Substring(1);
            }

            //---
            // Calculate variant type
            //---
            if (ReferenceAlleleString == AlternateAlleleString)
                this.VarType = VariantType.INTERVAL;
            else if (ReferenceAlleleString.Length == 1 && AlternateAlleleString.Length == 1)
                this.VarType = VariantType.SNV;
            else if (ReferenceAlleleString.Length == AlternateAlleleString.Length)
                this.VarType = VariantType.MNP;
            else if (ReferenceAlleleString.Length < AlternateAlleleString.Length && AlternateAlleleString.StartsWith(ReferenceAlleleString))
                this.VarType = VariantType.INS;
            else if (ReferenceAlleleString.Length > AlternateAlleleString.Length && ReferenceAlleleString.StartsWith(AlternateAlleleString))
                this.VarType = VariantType.DEL;
            else
                this.VarType = VariantType.MIXED;

            //---
            // Start and end position
            // 	- Start is always the leftmost base
            //	- End is always the rightmost affected base in the ReferenceAlleleStringerence genome
            //---
            if (VarType == VariantType.INS || VarType == VariantType.SNV)
            {
                // These changes only affect one position in the reference genome
                OneBasedEnd = OneBasedStart;
            }
            else // if (isDel() || isMnp()) {
            {
                // Update 'end' position
                if (ReferenceAlleleString.Length > 1) OneBasedEnd = OneBasedStart + ReferenceAlleleString.Length - 1;
            }
        }

        public bool isDel()
        {
            return (VarType == VariantType.DEL);
        }

        public bool isDup()
        {
            return (VarType == VariantType.DUP);
        }

        public bool isElongation()
        {
            return LengthChange() > 0;
        }

        public bool isInDel()
        {
            return (VarType == VariantType.INS) || (VarType == VariantType.DEL);
        }

        public bool isIns()
        {
            return VarType == VariantType.INS;
        }

        public bool isInterval()
        {
            return VarType == VariantType.INTERVAL;
        }

        public bool isInv()
        {
            return VarType == VariantType.INV;
        }

        public bool isMixed()
        {
            return VarType == VariantType.MIXED;
        }

        public bool isMnv()
        {
            return VarType == VariantType.MNP;
        }

        public bool isSnv()
        {
            return VarType == VariantType.SNV;
        }

        public bool isStructural()
        {
            return isDel() || isInv() || isDup();
        }

        /// <summary>
        /// Calculate the number of bases of change in length
        /// </summary>
        /// <returns></returns>
        public long LengthChange()
        {
            if (isSnv() || isMnv())
                return 0;

            // This is a length changing Variant (i.e. Insertions, deletion, or mixed change)
            // Calculate the number of bases of change in length
            if (ReferenceAlleleString != "" || AlternateAlleleString != "") return AlternateAlleleString.Length - ReferenceAlleleString.Length;

            // Default to traditional apporach for imprecise and structural variants
            return OneBasedEnd - OneBasedStart;
        }

        public string NetChange(bool reverseStrand)
        {
            if (isDel())
            {
                return reverseStrand ? SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, ReferenceAlleleString).GetReverseComplementedSequence()) : ReferenceAlleleString; // Deletion have empty 'alt'
            }
            return reverseStrand ? SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, AlternateAlleleString).GetReverseComplementedSequence()) : AlternateAlleleString;
        }

        public string NetChange(Interval interval)
        {
            string netChange = AlternateAlleleString;
            if (isDel()) netChange = ReferenceAlleleString; // In deletions 'alt' is empty

            long removeBefore = interval.OneBasedStart - OneBasedStart;
            if (removeBefore > 0)
            {
                if (removeBefore >= netChange.Length) return ""; // Nothing left
            }
            else removeBefore = 0;

            long removeAfter = OneBasedEnd - interval.OneBasedEnd;
            if (removeAfter > 0)
            {
                if ((removeBefore + removeAfter) >= netChange.Length) return ""; // Nothing left
            }
            else removeAfter = 0;

            // Remove leading and trailing parts
            netChange = netChange.Substring((int)removeBefore, netChange.Length - (int)removeAfter);

            return netChange;
        }

        #endregion Public Methods
    }
}