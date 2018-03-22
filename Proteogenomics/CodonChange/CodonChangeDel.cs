using Bio.Extensions;
using System;
using System.Collections.Generic;

namespace Proteogenomics
{
    public class CodonChangeDel
        : CodonChangeStructural
    {
        public CodonChangeDel(Variant variant, Transcript transcript, List<VariantEffect> variantEffects)
            : base(variant, transcript, variantEffects)
        {
            ReturnNow = false;
            RequireNetCdsChange = true;
        }

        /// <summary>
        /// Analyze deletions in this transcript.
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        protected override bool ChangeCodon(Exon exon)
        {
            // Is there any net effect?
            if (NetCodingSequenceChange == "") return false;

            EffectType effType = EffectType.NONE;

            if (Variant.Includes(exon))
            {
                /**
                 * An exon has been entirely removed
                 */
                CodonsReference = "";
                CodonsAlternate = "";
                CodonStartNumber = CodonStartIndex = -1;
                effType = EffectType.EXON_DELETED;
            }
            else if (NetCodingSequenceChange.Length % CODON_SIZE != 0)
            {
                /**
                 * Length not multiple of CODON_SIZE => FRAME_SHIFT
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Delete 'AA' pos 0:	ACC CGG GAA ACC CGG GAA ACC CGG G
                 * 		Delete 'AA' pos 1:	ACC CGG GAA ACC CGG GAA ACC CGG G
                 * 		Delete 'AC' pos 2:	AAC CGG GAA ACC CGG GAA ACC CGG G
                 */
                CodonsReference = CodonsRef();
                CodonsAlternate = "";
                effType = EffectType.FRAME_SHIFT;
            }
            else if (CodonStartIndex == 0)
            {
                /**
                 * Length multiple of CODON_SIZE and insertion happens at codon boundary => CODON_DELETION
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Delete 'AAA' pos 0:	CCC GGG AAA CCC GGG AAA CCC GGG
                 */
                CodonsReference = CodonsRef();
                CodonsAlternate = "";
                effType = EffectType.CODON_DELETION;
            }
            else
            {
                /**
                 * Length multiple of CODON_SIZE and insertion does not happen at codon boundary => CODON_CHANGE_PLUS_CODON_DELETION
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Delete 'AAC' pos 1:	ACC GGG AAA CCC GGG AAA CCC GGG
                 * 		Delete 'ACC' pos 2:	AAC GGG AAA CCC GGG AAA CCC GGG
                 */
                CodonsReference = CodonsRef();
                CodonsAlternate = CodonsAlt();

                if (CodonsAlternate == "" || CodonsReference.StartsWith(CodonsAlternate))
                {
                    /**
                     * Note: It might happen that the last codon of the exon was deleted.
                     *       In this case there is no 'CODON_CHANGE'
                     * E.g.
                     * 		Original:				AAA CCC GGG AAA CCC GGG AAA CCC GGG
                     * 		Delete 'GGG' pos 24:	ACC CCC GGG AAA CCC GGG AAA CCC
                     *
                     * Note2: It may also be the case that the deleted bases are equal to the following ones.
                     *  E.g.
                     *  	Original:			ACG TCG TCC GGG AAA CCC GGG AAA CCC GGG
                     *  	Delete 'CGT' pos 1:	ACG TCC GGG AAA CCC GGG AAA CCC GGG
                     */
                    effType = EffectType.CODON_DELETION;
                }
                else
                {
                    effType = EffectType.CODON_CHANGE_PLUS_CODON_DELETION;
                }
            }

            Effect(exon, effType, false);

            return true;
        }

        /// <summary>
        /// Get alternative codons
        /// </summary>
        /// <returns></returns>
        protected override string CodonsAlt()
        {
            if (NetCodingSequenceChange == "") return "";

            int after = NetCodingSequenceChange.Length + CodonStartIndex;

            string prefix = CodonsReference.Length >= CodonStartIndex ? CodonsReference.Substring(0, CodonStartIndex) : CodonsReference;
            string suffix = CodonsReference.Length > after ? CodonsReference.Substring(after) : "";

            string codonsAlt = prefix + suffix;
            return codonsAlt;
        }

        /// <summary>
        /// Get original codons in CDS
        /// </summary>
        /// <returns></returns>
        protected override string CodonsRef()
        {
            if (NetCodingSequenceChange == "") return "";

            long min = Variant.OneBasedStart;
            long max = Variant.OneBasedEnd;
            long cdsBaseMin = CdsBaseNumber(min);
            long cdsBaseMax = CdsBaseNumber(max);

            // Swap?
            if (Transcript.Strand != "+")
            {
                long swap = cdsBaseMin;
                cdsBaseMin = cdsBaseMax;
                cdsBaseMax = swap;
            }

            if (cdsBaseMax < cdsBaseMin) throw new Exception("This should never happen!\n\tcdsBaseMin: " + cdsBaseMin + "\n\tcdsBaseMax: " + cdsBaseMax + "\n\tmin: " + min + "\n\tmax: " + max + "\n\tSeqChange: " + Variant + "\n\ttranscript: " + Transcript + "\n\tCDS.len: " + Transcript.RetrieveCodingSequence().Count);

            long maxCodon = cdsBaseMax / CODON_SIZE;
            long minCodon = cdsBaseMin / CODON_SIZE;
            long oldCodonCdsStart = (CODON_SIZE * minCodon);
            long oldCodonCdsEnd = (CODON_SIZE * (maxCodon + 1)) - 1;

            string codons = "";
            if (oldCodonCdsEnd >= Transcript.RetrieveCodingSequence().Count) codons = SequenceExtensions.ConvertToString(Transcript.RetrieveCodingSequence()).Substring((int)oldCodonCdsStart);
            else codons = SequenceExtensions.ConvertToString(Transcript.RetrieveCodingSequence().GetSubSequence(oldCodonCdsStart, oldCodonCdsEnd + 1));

            return codons;
        }

        protected override void EffectTranscript()
        {
            EffectNoCodon(Transcript, EffectType.TRANSCRIPT_DELETED);
        }

        /// <summary>
        /// Whole exon/s deleted?
        /// </summary>
        private void ExonLoss()
        {
            foreach (Exon ex in Transcript.Exons)
            {
                if (Variant.Includes(ex)) EffectNoCodon(ex, EffectType.EXON_DELETED);
            }
        }

        /// <summary>
        /// Deletion analysis using full transcript information. This is done only when the variant affects more than one exons.
        /// </summary>
        protected override void Exons()
        {
            if (exonFull == 0 && exonPartial == 1)
            {
                // Variant partially affects only one exon?
                // => Use the standard (by exon) method
                codonChangeSuper();
                return;
            }
            else if (exonFull > 0)
            {
                // Full exons deleted
                ExonLoss();

                // Only whole exons deleted? We are done
                if (exonPartial == 0) { return; }
            }

            //---
            // A combination of partial and full exons affected
            //---
            CodonsRefAlt();
            EffectType effType = EffectType.NONE;
            int lenDiff = cdsAlt.Length - cdsRef.Length;
            if (lenDiff % CODON_SIZE != 0)
            {
                effType = EffectType.FRAME_SHIFT;
            }
            else if (CodonStartIndex == 0)
            {
                effType = EffectType.CODON_DELETION;
            }
            else
            {
                if (CodonsAlternate == "" || CodonsReference.StartsWith(CodonsAlternate))
                {
                    effType = EffectType.CODON_DELETION;
                }
                else
                {
                    effType = EffectType.CODON_CHANGE_PLUS_CODON_DELETION;
                }
            }

            // Assign to first exon
            foreach (Exon ex in Transcript.Exons)
            {
                if (Variant.Includes(ex) || Variant.Intersects(ex))
                {
                    Exon = ex;
                    break;
                }
            }

            // Add variant effect
            Effect(Exon, effType, false);
        }

        protected override void ExonsCoding()
        {
            // Nothing to do
        }

        protected override void ExonsNoncoding()
        {
            if (exonFull > 0) EffectNoCodon(Transcript, EffectType.EXON_DELETED, EffectImpact.MODIFIER);
            if (exonPartial > 0) EffectNoCodon(Transcript, EffectType.EXON_DELETED_PARTIAL, EffectImpact.MODIFIER);
        }

        protected override void Intron()
        {
            EffectNoCodon(Transcript, EffectType.INTRON);
        }
    }
}