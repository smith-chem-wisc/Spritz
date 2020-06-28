using Bio;
using Bio.Extensions;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class UTR5Prime
        : UTR
    {
        private List<UTR5Prime> UTRs { get; set; }

        public UTR5Prime(Exon parent, string chromID, string source, string strand, long oneBasedStart, long oneBasedEnd) :
            base(parent, chromID, source, strand, oneBasedStart, oneBasedEnd)
        {
        }

        public UTR5Prime(UTR5Prime utr) : base(utr)
        {
        }

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "five_prime_UTR";

        public override bool is3Prime()
        {
            return false;
        }

        public override bool is5Prime()
        {
            return true;
        }

        public List<UTR5Prime> get5primeUtrs()
        {
            if (UTRs == null)
            {
                Transcript tr = (Transcript)FindParent(typeof(Transcript));

                // Get UTRs and sort them
                UTRs = tr.UTRs.OfType<UTR5Prime>().ToList();
                if (IsStrandPlus())
                {
                    UTRs = UTRs.OrderBy(u => u.OneBasedStart).ToList(); // Sort by start position
                }
                else
                {
                    UTRs = UTRs.OrderByDescending(u => u.OneBasedEnd).ToList(); // Sort by end position (reversed)
                }
            }

            return UTRs;
        }

        public string getSequence()
        {
            // Create UTR sequence
            StringBuilder sb = new StringBuilder();
            foreach (UTR5Prime utr in get5primeUtrs())
            {
                Exon ex = (Exon)utr.Parent;
                ISequence utrSeq = ex.Sequence;
                if (utr.Length() < utrSeq.Count) { utrSeq = utrSeq.GetSubSequence(0, utr.Length()); } // UTR5' may stop before end of exon
                sb.Append(SequenceExtensions.ConvertToString(utrSeq));
            }
            return sb.ToString();
        }

        /// <summary>
        /// Is a new start codon produced? @return New start codon (or empty string if there is no new start codon)
        /// </summary>
        /// <param name="chars"></param>
        /// <param name="pos"></param>
        /// <returns></returns>
        private static string startGained(char[] chars, long pos)
        {
            // Analyze all frames
            for (long i = Math.Max(0, pos - 2); i <= pos && i + 2 < chars.Length; i++)
            {
                string codon = "" + chars[i] + chars[i + 1] + chars[i + 2];
                if (CodonsStandard.START_CODONS.Contains(codon.ToUpper(CultureInfo.InvariantCulture)))
                {
                    return codon.ToUpper(CultureInfo.InvariantCulture); // This frame has a start codon?
                }
            }
            return "";
        }
    }
}