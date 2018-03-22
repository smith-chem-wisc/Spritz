using System.Collections.Generic;

namespace Proteogenomics
{
    public class CodonChangeIns
        : CodonChange
    {
        public CodonChangeIns(Variant seqChange, Transcript transcript, List<VariantEffect> changeEffects)
            : base(seqChange, transcript, changeEffects)
        {
            ReturnNow = true; // An insertion can only affect one exon
        }

        /// <summary>
        /// Analyze insertions in this transcript. Add changeEffect to 'changeEffect'
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        protected override bool ChangeCodon(Exon exon)
        {
            string netChange = Variant.NetChange(Transcript.Strand != "+");

            CodonsReference = CodonsRef();
            CodonsAlternate = CodonsAlt();

            EffectType effType = EffectType.NONE;

            if (netChange.Length % CODON_SIZE != 0)
            {
                /**
                 * Length not multiple of CODON_SIZE => FRAME_SHIFT
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Insert 'TT' pos 0:	TTA AAC CCG GGA AAC CCG GGA AAC CCG GG
                 * 		Insert 'TT' pos 1:	ATT AAC CCG GGA AAC CCG GGA AAC CCG GG
                 * 		Insert 'TT' pos 2:	AAT TAC CCG GGA AAC CCG GGA AAC CCG GG
                 */
                effType = EffectType.FRAME_SHIFT;
            }
            else if (CodonStartIndex == 0)
            {
                /**
                 * Length multiple of CODON_SIZE and insertion happens at codon boundary => CODON_INSERTION
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Insert 'TTT' pos 0:	TTT AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 */
                effType = EffectType.CODON_INSERTION;
            }
            else
            {
                /**
                 * Length multiple of CODON_SIZE and insertion does not happen at codon boundary => CODON_CHANGE_PLUS_CODON_INSERTION
                 * 	E.g. :
                 * 		Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Insert 'TTT' pos 1:	ATT TAA CCC GGG AAA CCC GGG AAA CCC GGG
                 * 		Insert 'TTT' pos 2:	AAT TTA CCC GGG AAA CCC GGG AAA CCC GGG
                 */
                if (CodonsAlternate.ToUpper().StartsWith(CodonsReference.ToUpper()))
                {
                    /**
                     *  May be the inserted base are equal to the old ones.
                     *  E.g.
                     *  	Original:			AAA CCC GGG AAA CCC GGG AAA CCC GGG
                     *  	Insert 'AAA' pos 1:	AAA AAA CCC GGG AAA CCC GGG AAA CCC GGG
                     */
                    effType = EffectType.CODON_INSERTION;
                }
                else
                {
                    effType = EffectType.CODON_CHANGE_PLUS_CODON_INSERTION;
                }
            }

            Effect(exon, effType, false);

            return true;
        }

        /// <summary>
        /// Get new (modified) codons
        /// </summary>
        /// <returns></returns>
        protected override string CodonsAlt()
        {
            // Inserts BEFORE base:
            //		- In positive strand that is BEFORE pos
            //		- In negative strand, that is AFTER pos
            int idx = CodonStartIndex + (Transcript.Strand != "+" ? 1 : 0);

            // Insertion: Concatenate...
            string prefix = CodonsReference.Length >= idx ? CodonsReference.Substring(0, idx) : CodonsReference; // First part of the codon
            string netChange = Variant.NetChange(Transcript.Strand != "+"); // Insertion
            string suffix = CodonsReference.Length >= idx ? CodonsReference.Substring(idx) : ""; // last part of the codon

            // New codon
            string codonsNew = prefix + netChange + suffix;

            return codonsNew;
        }
    }
}