using Bio;
using Bio.Algorithms.Translation;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class Translation
    {

        #region One-Frame Translation

        public static Protein OneFrameTranslation(TranscriptPossiblyWithVariants transcript)
        {
            ISequence dnaSequence = transcript.Sequence;
            ISequence rnaSequence = Transcription.Transcribe(transcript.GetExonsUsedInDerivation()[0].Strand == "+" ? dnaSequence : dnaSequence.GetReverseComplementedSequence());
            ISequence proteinSequence = ProteinTranslation.Translate(rnaSequence);
            string proteinBases = SequenceExtensions.ConvertToString(proteinSequence).Split('*')[0];
            return new Protein(proteinBases, transcript.ProteinID, null, null, null, transcript.ProteinAnnotation, transcript.ProteinAnnotation);
        }

        #endregion One-Frame Translation

        #region Three-Frame Translation

        /// <summary>
        /// Not used or tested right now...
        /// </summary>
        /// <param name="exons"></param>
        /// <param name="proteinID"></param>
        /// <returns></returns>
        public static Protein ThreeFrameTranslation(List<Exon> exons, string proteinID)
        {
            string seq = String.Join("", exons.Select(x => SequenceExtensions.ConvertToString(x.Sequence)));
            if (seq.Contains('N')) return null;
            ISequence dna_seq = new Sequence(Alphabets.DNA, seq);
            ISequence rna_seq = Transcription.Transcribe(exons[0].Strand == "+" ? dna_seq : dna_seq.GetReverseComplementedSequence());
            ISequence[] prot_seq = Enumerable.Range(0, 3).Select(i => ProteinTranslation.Translate(rna_seq, i)).ToArray();

            //return the protein sequence corresponding to the longest ORF
            return new Protein(prot_seq.SelectMany(s => SequenceExtensions.ConvertToString(s).Split('*')).OrderByDescending(s => s.Length).FirstOrDefault(), proteinID);
        }

        #endregion Three-Frame Translation

    }
}
