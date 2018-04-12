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
            , MNV // Multiple nucleotide polymorphism (i.e. several bases are changed)
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
        /// Chromosme for this variant
        /// </summary>
        public Chromosome Chromosome { get; set; }

        /// <summary>
        /// Type of this variant
        /// </summary>
        public VariantType VarType { get; set; }

        /// <summary>
        /// The reference allele
        /// </summary>
        public Allele ReferenceAllele { get; set; }

        /// <summary>
        /// The reference allele as a string
        /// </summary>
        public string ReferenceAlleleString { get; set; }

        /// <summary>
        /// The first allele
        /// </summary>
        public Allele FirstAllele { get; set; }

        /// <summary>
        /// The alternate allele as a string
        /// </summary>
        public string FirstAlleleString { get; set; }

        /// <summary>
        /// The first allele
        /// </summary>
        public Allele SecondAllele { get; set; }

        /// <summary>
        /// The alternate allele as a string
        /// </summary>
        public string SecondAlleleString { get; set; }

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
        public int FirstAlleleDepth { get; set; }

        /// <summary>
        /// Depth of reads across this position containing the alternate allele
        /// </summary>
        public int SecondAlleleDepth { get; set; }

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

        public Variant(Interval parent, VariantContext variant, Genome genome)
            : this(parent, variant, genome.Chromosomes.FirstOrDefault(c => c.FriendlyName == variant.Chr))
        {
        }

        public Variant(Interval parent, VariantContext variant, Chromosome chromosome)
            : base(parent, variant.Chr, "+", variant.Start, variant.End, null)
        {
            Chromosome = chromosome;
            ReferenceAllele = variant.Reference;
            ReferenceAlleleString = variant.Reference.BaseString;
            Genotype = variant.Genotypes.First(); // assumes just one sample
            GenotypeType = Genotype.Type;
            ReadDepth = Genotype.DP;
            FirstAlleleDepth = Genotype.AD == null ? -1 : Genotype.AD[0];
            SecondAlleleDepth = Genotype.AD == null ? -1 : Genotype.AD[1];
            FirstAllele = Genotype.GetAllele(0);
            FirstAlleleString = FirstAllele.BaseString;
            SecondAllele = Genotype.GetAllele(1);
            SecondAlleleString = SecondAllele.BaseString;
            CalculateType();

            if (variant.Attributes.TryGetValue("ANN", out object snpEffAnnotationObj))
            {
                if (snpEffAnnotationObj as string[] != null) { SnpEffAnnotations = (snpEffAnnotationObj as string[]).Select(x => new SnpEffAnnotation(this, x)).ToList(); }
                if (snpEffAnnotationObj as string != null) { SnpEffAnnotations = new List<SnpEffAnnotation> { new SnpEffAnnotation(this, snpEffAnnotationObj as string) }; }
            }
        }

        #endregion Public Constructor

        #region Public Methods

        public void CalculateType()
        {
            // Sanity check
            if (SecondAlleleString.IndexOf(',') >= 0 || SecondAlleleString.IndexOf('/') >= 0)
            {
                throw new ArgumentException("Variants with multiple ALTs are not allowed (ALT: '" + SecondAlleleString + "')");
            }

            // Remove leading char (we still have some test cases using old TXT format)
            if (ReferenceAlleleString == "*") { ReferenceAlleleString = ""; }

            // Possibly some old formatting with +, -, and =
            if (SecondAlleleString.StartsWith("+"))
            {
                // Insertion
                SecondAlleleString = ReferenceAlleleString + SecondAlleleString.Substring(1);
            }
            else if (SecondAlleleString.StartsWith("-"))
            {
                // Deletion
                ReferenceAlleleString = SecondAlleleString.Substring(1);
                SecondAlleleString = "";
            }
            else if (SecondAlleleString.StartsWith("="))
            {
                // Mixed variant
                SecondAlleleString = SecondAlleleString.Substring(1);
            }

            //---
            // Calculate variant type
            //---
            if (ReferenceAlleleString == SecondAlleleString)
            {
                VarType = VariantType.INTERVAL;
            }
            else if (ReferenceAlleleString.Length == 1 && SecondAlleleString.Length == 1)
            {
                VarType = VariantType.SNV;
            }
            else if (ReferenceAlleleString.Length == SecondAlleleString.Length)
            {
                VarType = VariantType.MNV;
            }
            else if (ReferenceAlleleString.Length < SecondAlleleString.Length && SecondAlleleString.StartsWith(ReferenceAlleleString))
            {
                VarType = VariantType.INS;
            }
            else if (ReferenceAlleleString.Length > SecondAlleleString.Length && ReferenceAlleleString.StartsWith(SecondAlleleString))
            {
                VarType = VariantType.DEL;
            }
            else
            {
                VarType = VariantType.MIXED;
            }

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
            return VarType == VariantType.DEL;
        }

        public bool isDup()
        {
            return VarType == VariantType.DUP;
        }

        public bool isElongation()
        {
            return LengthChange() > 0;
        }

        public bool isInDel()
        {
            return VarType == VariantType.INS || VarType == VariantType.DEL;
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
            return VarType == VariantType.MNV;
        }

        public bool isSnv()
        {
            return VarType == VariantType.SNV;
        }

        public bool isStructural()
        {
            return isDel() || isInv() || isDup();
        }

        public bool isVariant()
        {
            return VarType != VariantType.INTERVAL;
        }

        /// <summary>
        /// Calculate the number of bases of change in length
        /// </summary>
        /// <returns></returns>
        public long LengthChange()
        {
            if (isSnv() || isMnv())
            {
                return 0;
            }

            // This is a length changing Variant (i.e. Insertions, deletion, or mixed change)
            // Calculate the number of bases of change in length
            if (ReferenceAlleleString != "" || SecondAlleleString != "") { return SecondAlleleString.Length - ReferenceAlleleString.Length; }

            // Default to traditional apporach for imprecise and structural variants
            return OneBasedEnd - OneBasedStart;
        }

        public string NetChange(bool reverseStrand)
        {
            if (isDel())
            {
                return reverseStrand ? 
                    SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, ReferenceAlleleString).GetReverseComplementedSequence()) : 
                    ReferenceAlleleString; // Deletion have empty 'alt'
            }
            return reverseStrand ?
                SequenceExtensions.ConvertToString(new Sequence(Alphabets.DNA, SecondAlleleString).GetReverseComplementedSequence()) : 
                SecondAlleleString;
        }

        public string NetChange(Interval interval)
        {
            string netChange = isDel() ? ReferenceAlleleString : SecondAlleleString; // In deletions 'alt' is empty

            long removeBefore = interval.OneBasedStart - OneBasedStart;
            if (removeBefore > 0)
            {
                if (removeBefore >= netChange.Length) { return ""; }// Nothing left
            }
            else
            {
                removeBefore = 0;
            }

            long removeAfter = OneBasedEnd - interval.OneBasedEnd;
            if (removeAfter > 0)
            {
                if ((removeBefore + removeAfter) >= netChange.Length) { return ""; } // Nothing left
            }
            else
            {
                removeAfter = 0;
            }

            // Remove leading and trailing parts
            int netChangeLen = netChange.Length - (int)removeAfter - (int)removeBefore;
            netChange = netChange.Substring((int)removeBefore, netChangeLen);

            return netChange;
        }

        public override int CompareTo(Interval i2)
        {
            int comp = base.CompareTo(i2);
            if (comp != 0) { return comp; }

            Variant v2 = i2 as Variant;
            if (v2 == null) { return GetType().Name.CompareTo(i2.GetType().Name); }

            comp = ReferenceAlleleString.CompareTo(v2.ReferenceAlleleString);
            if (comp != 0) { return comp; }

            comp = SecondAlleleString.CompareTo(v2.SecondAlleleString);
            if (comp != 0) { return comp; }

            return 0;
        }

        public override string ToString()
        {
            if ((ReferenceAlleleString == null || ReferenceAlleleString == "") && (SecondAlleleString == null || SecondAlleleString == ""))
            {
                return ChromosomeID + ":" + OneBasedStart.ToString() + "-" + OneBasedEnd.ToString() + "[" + VarType.ToString() + "]";
            }

            return ChromosomeID //
                    + ":" + OneBasedStart.ToString() //
                    + "_" + ReferenceAlleleString //
                    + "/" + SecondAlleleString
                    + "_" + Genotype.Type.ToString();
        }

        #endregion Public Methods
    }
}