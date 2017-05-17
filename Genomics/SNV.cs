using System.Collections.Generic;

namespace Genomics
{
    public class SNV : SequenceVariant
    {

        #region Public Constructor

        public SNV(Chromosome chrom, int position, string id, string reference, string alternate, double qual, string filter, Dictionary<string, string> info)
            : base(chrom, position, id, reference, alternate, qual, filter, info)
        { }

        #endregion Public Constructor

        #region Public Methods

        public bool is_missense()
        {
            return false;
        }

        public bool is_start_gain()
        {
            return false;
        }

        public bool is_stop_loss()
        {
            return false;
        }

        public bool parse_aa_change(string aa_change)
        {
            //aa_abbrev_dict = amino_acids_3to1
            //aa_change_regex = '([A-Z])(\d+)([A-Z])'  # G528R
            //aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))'  # p.Gly528Arg/c.1582G>C
            //aa_pos = None  # 1-based position
            //ref_aa, alt_aa = '_', '_'
            //m = re.match(aa_change_regex, aa_change)  # parse aa_change, and get AA change position and alternate Animo Acid
            //if m:
            //    aa_pos = int(m.groups()[1])
            //    ref_aa = m.groups()[0]
            //    alt_aa = m.groups()[2]
            //else:
            //    m = re.match(aa_hgvs_regex, aa_change)
            //    if m:
            //        aa_pos = int(m.groups()[1])
            //        ref_aa = aa_abbrev_dict[m.groups()[0]]
            //        alt_aa = aa_abbrev_dict[m.groups()[2]]
            //return aa_pos, ref_aa, alt_aa
            return false;
        }

        #endregion Public Methods

    }
}
