using Bio;
using Bio.Extensions;
using System;
using System.Globalization;

namespace Proteogenomics
{
    public class CodonChangeSnv
        : CodonChange
    {
        public CodonChangeSnv(Variant variant, Transcript transcript, VariantEffects variantEffects)
            : base(variant, transcript, variantEffects)
        {
            ReturnNow = true; // A SNP can only affect one exon
        }

        /// <summary>
        /// Analyze SNPs in this transcript. Add changeEffect to 'changeEffect'
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        protected override bool ChangeCodon(Exon exon)
        {
            // Get old and new codons
            CodonsReference = CodonsRef();
            CodonsAlternate = CodonsAlt();

            // Use a generic low priority variant, this allows 'AdditionalEffects' to override it
            Effect(exon, EffectType.CODON_CHANGE, true);

            if (CodonsReference == "")
            {
                VariantEffects.AddErrorWarning(Variant, ErrorWarningType.ERROR_MISSING_CDS_SEQUENCE);
            }

            return true;
        }

        /// <summary>
        /// Get new (modified) codons
        /// </summary>
        /// <returns></returns>
        protected override string CodonsAlt()
        {
            // Was there a problem getting 'codonsOld'? => We cannot do anything
            if (CodonsReference == "")
            {
                return "";
            }

            char[] codonChars = CodonsReference.ToLower(CultureInfo.InvariantCulture).ToCharArray();
            char snpBase = Variant.NetChange(Transcript.IsStrandMinus())[0];
            if (CodonStartIndex < codonChars.Length)
            {
                codonChars[CodonStartIndex] = char.ToUpper(snpBase);
            }

            string codonsNew = new string(codonChars);
            return codonsNew;
        }

        /// <summary>
        /// Get original codons in CDS
        /// </summary>
        /// <returns></returns>
        protected override string CodonsRef()
        {
            int numCodons = 1;

            // Get CDS
            ISequence cdsStr = Transcript.RetrieveCodingSequence();
            long cdsLen = cdsStr.Count;

            // Calculate minBase (first codon base in the CDS)
            int minBase = CodonStartNumber * CODON_SIZE;
            if (minBase < 0) { minBase = 0; }

            // Calculate maxBase (last codon base in the CDS)
            long maxBase = CodonStartNumber * CODON_SIZE + numCodons * CODON_SIZE;
            if (maxBase > cdsLen) { maxBase = cdsLen; }

            // Sanity checks
            if (cdsLen == 0 // Empty CDS => Cannot get codon (e.g. one or more exons are missing their sequences
                    || (cdsLen <= minBase) // Codon past CDS sequence => Cannot get codon
            ) return "";

            // Create codon sequence
            char[] codonChars = SequenceExtensions.ConvertToString(cdsStr).Substring(minBase, CODON_SIZE).ToLower(CultureInfo.InvariantCulture).ToCharArray();

            // Capitatlize changed base
            if (CodonStartIndex < codonChars.Length) { codonChars[CodonStartIndex] = char.ToUpper(codonChars[CodonStartIndex]); }
            string codon = new String(codonChars);

            return codon;
        }
    }
}