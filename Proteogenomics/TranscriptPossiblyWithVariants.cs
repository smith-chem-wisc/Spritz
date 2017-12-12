using System.Collections.Generic;
using System.Linq;
using Bio;

namespace Proteogenomics
{
    public class TranscriptPossiblyWithVariants
    {

        #region Public Properties

        public Transcript Transcript { get; set; }

        public ISequence Sequence { get; set; }

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
            Sequence = sequence;
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
            return Sequence.Count >= 3 && !ContainsAmbiguity;
        }

        #endregion Public Method

    }
}
