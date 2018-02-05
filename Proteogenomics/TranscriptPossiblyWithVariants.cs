using System.Collections.Generic;
using System.Linq;
using Bio;

namespace Proteogenomics
{
    public class TranscriptPossiblyWithVariants
    {

        #region Public Properties

        public Transcript Transcript { get; set; }

        public ISequence ReferenceTranscriptSequence { get; set; }

        public ISequence VariantTranscriptSequence { get; set; }

        public int[] Indices { get; set; }

        public bool Mitochondrial { get; set; }

        public bool ContainsAmbiguity { get; set; }

        public bool DerivedFromCodingDomainSequences { get; set; }

        public long ZeroBasedCodingStart { get; set; }
        
        public List<Variant> Variants { get; set; }

        public string ProteinID { get; set; }

        public string ProteinAnnotation { get; set; } // Used as the full name of the protein

        #endregion Public Properties

        #region Public Constructor

        public TranscriptPossiblyWithVariants(Transcript transcript, bool usedCodingDomainSequences, ISequence sequence, bool containsAmbiguity, List<Variant> variants)
        {
            Transcript = transcript;
            VariantTranscriptSequence = sequence;
            ContainsAmbiguity = containsAmbiguity;
            ZeroBasedCodingStart = 0;
            Variants = variants.OrderBy(v => v.OneBasedStart).ToList();
            ProteinID = transcript.ProteinID;
            DerivedFromCodingDomainSequences = usedCodingDomainSequences;
        }

        #endregion Public Constructor

        #region Public Method

        public List<Exon> GetExonsUsedInDerivation()
        {
            return DerivedFromCodingDomainSequences ? Transcript.CodingDomainSequences : Transcript.Exons;
        }

        public bool OkayToTranslate()
        {
            return VariantTranscriptSequence.Count >= 3 && !ContainsAmbiguity;
        }

        public void PrepareForTranslation()
        {
            List<Exon> originalExons = GetExonsUsedInDerivation();
            if (Transcript.Strand != "+") originalExons.Reverse();
            IAlphabet originalAlphabet = originalExons.OrderByDescending(x => x.Sequence.Alphabet.Count).First().Sequence.Alphabet;
            ReferenceTranscriptSequence = new Sequence(originalAlphabet, originalExons.SelectMany(x => Transcript.Strand == "+" ? x.Sequence : x.Sequence.GetReverseComplementedSequence()).ToArray());
            Indices = GetExonsUsedInDerivation().SelectMany(x =>
                Transcript.Strand == "+" ?
                    Enumerable.Range((int)x.OneBasedStart, (int)(x.OneBasedEnd - x.OneBasedStart + 1)) :
                    Enumerable.Range((int)x.OneBasedStart, (int)(x.OneBasedEnd - x.OneBasedStart + 1)).Reverse()
                ).ToArray();
            Mitochondrial = Transcript.Gene.ChromosomeID.StartsWith("M");
        }

        #endregion Public Method

    }
}
