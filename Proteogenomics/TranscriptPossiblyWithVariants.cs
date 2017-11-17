using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public class TranscriptPossiblyWithVariants
    {

        #region Public Properties

        public Transcript Transcript { get; set; }

        public string Sequence { get; set; }

        public bool DerivedFromCodingDomainSequences { get; set; }

        public long ZeroBasedCodingStart { get; set; }
        
        public List<Variant> Variants { get; set; }

        public string ProteinID { get; set; }

        public string ProteinAnnotation { get; set; } // Used as the full name of the protein

        #endregion Public Properties

        #region Public Constructor

        public TranscriptPossiblyWithVariants(Transcript transcript, bool usedCodingDomainSequences, string sequence, List<Variant> variants)
        {
            Transcript = transcript;
            Sequence = sequence;
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
            return Sequence.Length >= 3 && !Sequence.Contains('N');
        }

        #endregion Public Method

    }
}
