using Bio;
using Bio.Algorithms.Translation;
using Bio.Extensions;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteogenomics
{
    public class Transcript
        : Interval
    {
        /// <summary>
        /// Used to construct upstream and downstream reginos
        /// </summary>
        public static readonly int DEFAULT_UP_DOWN_LENGTH = 5000;

        /// <summary>
        /// Keeping track of transcript IDs leading to failures in executing combinitorics
        /// </summary>
        public static List<string> combinatoricFailures = new List<string>();

        private List<Exon> _ExonSortedStrand;
        private List<CDS> _CdsSortedStrand;
        private Protein protein;
        private ISequence _CodingSequence;
        private Exon FirstCodingExon;
        private Exon LastCodingExon;

        /// <summary>
        /// Constructor from the GFF3 reader information, including IDs, strand and Protein ID if available.
        /// </summary>
        /// <param name="id"></param>
        /// <param name="gene"></param>
        /// <param name="metadata"></param>
        /// <param name="ProteinID"></param>
        public Transcript(string id, Gene gene, string source, string strand, long oneBasedStart, long oneBasedEnd, string proteinID, MetadataListItem<List<string>> featureMetadata)
            : base(gene, gene.ChromosomeID, source, strand, oneBasedStart, oneBasedEnd)
        {
            ID = id;
            ProteinID = proteinID ?? id;
            Gene = gene;
            FeatureMetadata = featureMetadata;
        }

        /// <summary>
        /// Copy this transcript
        /// </summary>
        /// <param name="transcript"></param>
        public Transcript(Transcript transcript)
            : this(transcript.ID, transcript.Gene, transcript.Source, transcript.Strand,
                  transcript.OneBasedStart, transcript.OneBasedEnd, transcript.ProteinID, transcript.FeatureMetadata)
        {
            VariantAnnotations = new List<string>(transcript.VariantAnnotations);
            ProteinSequenceVariations = new HashSet<SequenceVariation>(transcript.ProteinSequenceVariations);
            Exons = new List<Exon>(transcript.Exons.Select(x => new Exon(this, x.Sequence, x.Source, x.OneBasedStart, x.OneBasedEnd, x.ChromosomeID, x.Strand, x.FeatureMetadata)));
            CodingDomainSequences = new List<CDS>(transcript.CodingDomainSequences.Select(cds => new CDS(this, cds.ChromosomeID, cds.Source, cds.Strand, cds.OneBasedStart, cds.OneBasedEnd, cds.StartFrame)));
            SetRegions(this);
        }

        /// <summary>
        /// The transcript ID
        /// </summary>
        public string ID { get; set; }

        /// <summary>
        /// The parent gene containing this transcript.
        /// </summary>
        public Gene Gene { get; set; }

        /// <summary>
        /// A list of exons, possibly including untranslated regions (UTRs).
        /// </summary>
        public List<Exon> Exons { get; set; } = new List<Exon>();

        /// <summary>
        /// Exons sorted by start position or by reverse end position if on reverse strand
        /// </summary>
        public List<Exon> ExonsSortedStrand
        {
            get
            {
                _ExonSortedStrand = _ExonSortedStrand ?? SortedStrand(Exons.OfType<Interval>().ToList()).OfType<Exon>().ToList();
                return _ExonSortedStrand;
            }
            private set
            {
                _ExonSortedStrand = value;
            }
        }

        /// <summary>
        /// Exons sorted by start position or by reverse end position if on reverse strand
        /// </summary>
        public List<CDS> CdsSortedStrand
        {
            get
            {
                _CdsSortedStrand = _CdsSortedStrand ?? SortedStrand(CodingDomainSequences.OfType<Interval>().ToList()).OfType<CDS>().ToList();
                return _CdsSortedStrand;
            }
            private set
            {
                _CdsSortedStrand = value;
            }
        }

        /// <summary>
        /// List of the untranslated regions
        /// </summary>
        public List<UTR> UTRs { get; set; } = new List<UTR>();

        /// <summary>
        /// Downstream region, 5 kb by default
        /// </summary>
        public Downstream Downstream { get; set; }

        /// <summary>
        /// Upstream region, 5 kb by default
        /// </summary>
        public Upstream Upstream { get; set; }

        /// <summary>
        /// Introns for this transcript
        /// </summary>
        public List<Intron> Introns { get; set; } = new List<Intron>();

        /// <summary>
        /// Coding domain sequence (CDS) information
        /// </summary>
        public List<CDS> CodingDomainSequences { get; set; } = new List<CDS>();

        /// <summary>
        /// Start of coding sequence
        /// </summary>
        public long CdsOneBasedStart { get; set; } = -1;

        /// <summary>
        /// End of coding sequence
        /// </summary>
        public long CdsOneBasedEnd { get; set; } = -1;

        public long[] Cds2Pos { get; set; }
        public long[] AA2Pos { get; set; }

        /// <summary>
        /// The protein ID derived from coding transcripts; this is imported from the GFF3 file if available.
        /// </summary>
        public string ProteinID { get; set; }

        /// <summary>
        /// Annotations for the variants applied to this transcript
        /// </summary>
        public List<string> VariantAnnotations { get; set; } = new List<string>();

        /// <summary>
        /// List of protein sequence variations for annotating XML files
        /// </summary>
        public HashSet<SequenceVariation> ProteinSequenceVariations { get; set; } = new HashSet<SequenceVariation>();

        /// <summary>
        /// Feature name used for writing GTF files
        /// </summary>
        public override string FeatureType { get; } = "transcript";

        public MetadataListItem<List<string>> FeatureMetadata { get; private set; }

        /// <summary>
        /// Sets relevant regions for a transcript
        /// </summary>
        /// <param name="transcript"></param>
        public static void SetRegions(Transcript transcript)
        {
            transcript.Introns = transcript.CreateIntrons();
            transcript.UTRs = transcript.CreateUTRs();
            var updown = transcript.CreateUpDown(transcript.Gene.Chromosome);
            transcript.Upstream = updown.OfType<Upstream>().FirstOrDefault();
            transcript.Downstream = updown.OfType<Downstream>().FirstOrDefault();
        }

        /// <summary>
        /// Find base at genomic coordinate 'pos'
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public ISequence BaseAt(int pos)
        {
            CalcCdsStartEnd();
            Exon ex = FindExon(pos);
            if (ex == null) { return null; }
            return ex.basesAt(pos - ex.OneBasedStart, 1);
        }

        /// <summary>
        /// Calculate distance from transcript start to a position
        /// mRNA is roughly the same than cDNA. Strictly speaking mRNA
        /// has a poly-A tail and 5'cap.
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public long BaseNumber2MRnaPos(long pos)
        {
            long count = 0;
            foreach (Exon x in ExonsSortedStrand)
            {
                if (x.Intersects(pos))
                {
                    // Intersect this exon? Calculate the number of bases from the beginning
                    long dist = IsStrandPlus() ?
                        pos - x.OneBasedStart :
                        x.OneBasedEnd - pos;

                    // Sanity check
                    if (dist < 0)
                    {
                        throw new ArgumentException("Negative distance for position " + pos + ". This should never happen!\n" + this);
                    }

                    return count + dist;
                }

                count += x.Length();
            }
            return -1;
        }

        /// <summary>
        /// Calculate base number in a CDS where 'pos' maps
        ///
        /// usePrevBaseIntron: When 'pos' is intronic this method returns:
        /// 			- if(usePrevBaseIntron== false)  => The first base in the exon after 'pos' (i.e.first coding base after intron)
        /// 			- if(usePrevBaseIntron== true)   => The last base in the exon before 'pos'  (i.e.last coding base before intron)
        ///
        /// </summary>
        /// <param name="pos"></param>
        /// <param name="usePrevBaseIntron"></param>
        /// <returns>Base number or '-1' if it does not map to a coding base</returns>
        public long BaseNumberCds(long pos, bool usePrevBaseIntron)
        {
            // Doesn't hit this transcript?
            if (!Intersects(pos)) { return -1; }

            // Is it in UTR instead of CDS?
            if (UTRs.Any(utr => utr.Intersects(pos))) { return -1; }

            // Calculate cdsStart and cdsEnd (if not already done)
            CalcCdsStartEnd();

            // All exons..
            long firstCdsBaseInExon = 0; // Where the exon maps to the CDS (i.e. which CDS base number does the first base in this exon maps to).
            foreach (Exon eint in ExonsSortedStrand)
            {
                if (eint.Intersects(pos))
                {
                    long cdsBaseInExon = IsStrandPlus() ? // cdsBaseInExon: base number relative to the beginning of the coding part of this exon (i.e. excluding 5'UTRs)
                        pos - Math.Max(eint.OneBasedStart, CdsOneBasedStart) :
                        Math.Min(eint.OneBasedEnd, CdsOneBasedStart) - pos;

                    cdsBaseInExon = Math.Max(0, cdsBaseInExon);

                    return firstCdsBaseInExon + cdsBaseInExon;
                }
                else
                {
                    // Before exon begins?
                    if (IsStrandPlus() && pos < eint.OneBasedStart // Before exon begins (positive strand)?
                            || IsStrandMinus() && pos > eint.OneBasedEnd) // Before exon begins (negative strand)?
                        return firstCdsBaseInExon - (usePrevBaseIntron ? 1 : 0);
                }

                firstCdsBaseInExon += IsStrandPlus() ?
                    Math.Max(0, eint.OneBasedEnd - Math.Max(eint.OneBasedStart, CdsOneBasedStart) + 1) :
                    Math.Max(0, Math.Min(CdsOneBasedStart, eint.OneBasedEnd) - eint.OneBasedStart + 1);
            }

            return firstCdsBaseInExon - 1;
        }

        /// <summary>
        /// Calculate chromosome position as function of CDS number
        /// </summary>
        /// <returns>An array mapping 'cds2pos[cdsBaseNumber] = chromosmalPos'</returns>
        public long[] BaseNumberCds2Pos()
        {
            if (Cds2Pos != null) { return Cds2Pos; }

            CalcCdsStartEnd();

            Cds2Pos = new long[RetrieveCodingSequence().Count];
            for (int i = 0; i < Cds2Pos.Length; i++)
            {
                Cds2Pos[i] = -1;
            }

            long cdsMin = Math.Min(CdsOneBasedStart, CdsOneBasedEnd);
            long cdsMax = Math.Max(CdsOneBasedStart, CdsOneBasedEnd);

            // For each exon, add CDS position to array
            int cdsBaseNum = 0;
            foreach (Exon exon in ExonsSortedStrand)
            {
                long min = IsStrandPlus() ? exon.OneBasedStart : exon.OneBasedEnd;
                int step = IsStrandPlus() ? 1 : -1;
                for (long pos = min; exon.Intersects(pos) && cdsBaseNum < Cds2Pos.Length; pos += step)
                {
                    if (cdsMin <= pos && pos <= cdsMax)
                    {
                        Cds2Pos[cdsBaseNum++] = pos;
                    }
                }
            }

            return Cds2Pos;
        }

        /// <summary>
        /// Calculate CDS start and CDS end
        /// </summary>
        private void CalcCdsStartEnd()
        {
            // Do we need to calculate these values?
            // Note: In circular genomes, one of cdsStart / cdsEnd might be less
            //       than zero (we must check both)
            if (CdsOneBasedStart < 0 && CdsOneBasedEnd < 0)
            {
                // Calculate coding start (after 5 prime UTR)

                if (UTRs.Count == 0)
                {
                    // No UTRs => Use all exons
                    CdsOneBasedStart = IsStrandPlus() ? OneBasedEnd : OneBasedStart; // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
                    CdsOneBasedEnd = IsStrandPlus() ? OneBasedStart : OneBasedEnd; // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)

                    foreach (Exon ex in Exons)
                    {
                        if (IsStrandPlus())
                        {
                            CdsOneBasedStart = Math.Min(CdsOneBasedStart, ex.OneBasedStart);
                            CdsOneBasedEnd = Math.Max(CdsOneBasedEnd, ex.OneBasedEnd);
                        }
                        else
                        {
                            CdsOneBasedStart = Math.Max(CdsOneBasedStart, ex.OneBasedEnd);
                            CdsOneBasedEnd = Math.Min(CdsOneBasedEnd, ex.OneBasedStart);
                        }
                    }
                }
                else
                {
                    // We have to take into account UTRs
                    CdsOneBasedStart = IsStrandPlus() ? OneBasedStart : OneBasedEnd; // cdsStart is the position of the first base in the CDS (i.e. the first base after all 5'UTR)
                    CdsOneBasedEnd = IsStrandPlus() ? OneBasedEnd : OneBasedStart; // cdsEnd is the position of the last base in the CDS (i.e. the first base before all 3'UTR)
                    long cdsStartNotExon = CdsOneBasedStart;

                    foreach (UTR utr in UTRs)
                    {
                        if (utr is UTR5Prime)
                        {
                            if (IsStrandPlus()) { CdsOneBasedStart = Math.Max(CdsOneBasedStart, utr.OneBasedEnd + 1); }
                            else { CdsOneBasedStart = Math.Min(CdsOneBasedStart, utr.OneBasedStart - 1); }
                        }
                        else if (utr is UTR3Prime)
                        {
                            if (IsStrandPlus()) { CdsOneBasedEnd = Math.Min(CdsOneBasedEnd, utr.OneBasedStart - 1); }
                            else { CdsOneBasedEnd = Math.Max(CdsOneBasedEnd, utr.OneBasedEnd + 1); }
                        }
                    }

                    // Make sure cdsStart and cdsEnd lie within an exon
                    if (IsStrandPlus())
                    {
                        CdsOneBasedStart = FirstExonPositionAfter(CdsOneBasedStart);
                        CdsOneBasedEnd = LastExonPositionBefore(CdsOneBasedEnd);
                    }
                    else
                    {
                        CdsOneBasedStart = LastExonPositionBefore(CdsOneBasedStart);
                        CdsOneBasedEnd = FirstExonPositionAfter(CdsOneBasedEnd);
                    }

                    // We were not able to find cdsStart & cdsEnd within exon limits.
                    // Probably there is something wrong with the database and the transcript does
                    // not have a single coding base (e.g. all of it is UTR).
                    if (CdsOneBasedStart < 0 || CdsOneBasedEnd < 0)
                    {
                        CdsOneBasedStart = CdsOneBasedEnd = cdsStartNotExon;
                    }
                }
            }
        }

        /// <summary>
        /// Create a marker of the coding region in this transcript
        /// </summary>
        /// <returns></returns>
        public Interval CdsMarker()
        {
            Interval interval = new Interval(this);
            interval.OneBasedStart = IsStrandPlus() ? CdsOneBasedStart : CdsOneBasedEnd;
            interval.OneBasedEnd = IsStrandPlus() ? CdsOneBasedEnd : CdsOneBasedStart;
            return interval;
        }

        /// <summary>
        /// Find a CDS that matches exactly the exon
        /// </summary>
        /// <param name="exon"></param>
        /// <returns></returns>
        public CDS FindCds(Exon exon)
        {
            return CodingDomainSequences.FirstOrDefault(c => exon.Includes(c));
        }

        /// <summary>
        /// Return the an exon that intersects 'pos'
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        public Exon FindExon(int pos)
        {
            return Exons.FirstOrDefault(x => x.Includes(pos));
        }

        /// <summary>
        /// Return an exon intersecting 'marker' (first exon found)
        /// </summary>
        /// <param name="marker"></param>
        /// <returns></returns>
        public Exon FindExon(Interval marker)
        {
            return Exons.FirstOrDefault(x => x.Intersects(marker));
        }

        /// <summary>
        /// Find the first position after 'pos' within an exon
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        private long FirstExonPositionAfter(long pos)
        {
            foreach (Exon ex in Exons.OrderBy(x => x.OneBasedStart))
            {
                if (pos <= ex.OneBasedStart) { return ex.OneBasedStart; }
                if (pos <= ex.OneBasedEnd) { return pos; }
            }

            Console.WriteLine("WARNING: Cannot find first exonic position after " + pos + " for transcript '" + ID + "'");
            return -1;
        }

        /// <summary>
        /// Find the last position before 'pos' within an exon
        /// </summary>
        /// <param name="pos"></param>
        /// <returns></returns>
        private long LastExonPositionBefore(long pos)
        {
            long last = -1;
            foreach (Exon ex in Exons.OrderBy(x => x.OneBasedStart))
            {
                if (pos < ex.OneBasedStart)
                {
                    // Nothing found?
                    if (last < 0)
                    {
                        Console.WriteLine("WARNING: Cannot find last exonic position before " + pos + " for transcript '" + ID + "'");
                        return -1;
                    }
                    return last;
                }
                else if (pos <= ex.OneBasedEnd)
                {
                    return pos;
                }
                last = ex.OneBasedEnd;
            }

            if (last < 0)
            {
                Console.WriteLine("WARNING: Cannot find last exonic position before " + pos + " for transcript '" + ID + "'");
            }
            return pos;
        }

        #region Translate from Variant Nucleotide Sequences Methods

        /// <summary>
        /// Get the coding sequence for this transcript.
        /// SnpEff keeps track of the UTRs to figure this out. I suppose that will work, now that I'm using the interval tree to dive down to change those ranges.
        /// </summary>
        /// <returns></returns>
        public ISequence RetrieveCodingSequence()
        {
            if (_CodingSequence != null)
            {
                return _CodingSequence;
            }

            // Concatenate all exons
            List<Exon> exons = ExonsSortedStrand;
            StringBuilder sequence = new StringBuilder();
            int utr5len = 0;
            int utr3len = 0;

            // 5 prime UTR length
            foreach (UTR utr in UTRs.OfType<UTR5Prime>())
            {
                utr5len += (int)utr.Length();
            }

            // Append all exon sequences
            IAlphabet alphabet = Alphabets.AmbiguousDNA;
            bool missingSequence = false;
            foreach (Exon exon in exons)
            {
                missingSequence |= exon.Sequence == null; // If there is no sequence, we are in trouble
                sequence.Append(SequenceExtensions.ConvertToString(exon.Sequence)); // reverse complemented for reverse strand during loading
            }

            if (missingSequence)
            {
                _CodingSequence = new Sequence(Alphabets.DNA, ""); // One or more exons does not have sequence. Nothing to do
            }
            else
            {
                // OK, all exons have sequences

                // 3 prime UTR length
                foreach (UTR utr in UTRs.OfType<UTR3Prime>())
                {
                    utr3len += (int)utr.Length();
                }

                // Cut 5 prime UTR and 3 prime UTR points
                string dnaSequence = sequence.ToString();
                int subEnd = dnaSequence.Length - utr3len;
                int subLen = subEnd - utr5len;

                if (utr5len > subEnd)
                {
                    _CodingSequence = new Sequence(Alphabets.DNA, "");
                }
                else
                {
                    _CodingSequence = new Sequence(alphabet, dnaSequence.Substring(utr5len, subLen));
                }
            }
            return _CodingSequence;
        }

        /// <summary>
        /// Get strands sorted by start (forward) or end and in reverse (reverse strand)
        /// </summary>
        /// <returns></returns>
        private List<Interval> SortedStrand(List<Interval> intervals)
        {
            return IsStrandPlus() ?
                intervals.OrderBy(x => x.OneBasedStart).ToList() :
                intervals.OrderByDescending(x => x.OneBasedEnd).ToList();
        }

        public ISequence SplicedRNA()
        {
            bool ambiguity = ExonsSortedStrand.Any(x => x.Sequence.Alphabet.HasAmbiguity);
            return new Sequence(ambiguity ? Alphabets.AmbiguousDNA : Alphabets.DNA, String.Join("", ExonsSortedStrand.Select(x => SequenceExtensions.ConvertToString(x.Sequence))));
        }

        public Protein Protein()
        {
            return protein ?? Protein(null);
        }

        public Protein Protein(Dictionary<string, string> selenocysteineContaining)
        {
            return protein ?? Protein(RetrieveCodingSequence(), selenocysteineContaining);
        }

        private Protein Protein(ISequence dnaSeq, Dictionary<string, string> selenocysteineContaining)
        {
            selenocysteineContaining = selenocysteineContaining != null ? selenocysteineContaining : new Dictionary<string, string>();
            bool hasSelenocysteine = selenocysteineContaining.TryGetValue(ProteinID, out string selenocysteineContainingSeq);
            HashSet<int> uIndices = !hasSelenocysteine ?
                new HashSet<int>() :
                new HashSet<int>(Enumerable.Range(0, selenocysteineContainingSeq.Length).Where(i => selenocysteineContainingSeq[i] == 'U'));

            // Translate protein sequence, and replace amber stop codons with selenocysteines where appropriate
            ISequence proteinSequence = Translation.OneFrameTranslation(dnaSeq, Gene.Chromosome.Mitochondrial);
            string proteinBases = !hasSelenocysteine ?
                SequenceExtensions.ConvertToString(proteinSequence) :
                new string(Enumerable.Range(0, (int)proteinSequence.Count).Select(i => uIndices.Contains(i) && proteinSequence[i] == Alphabets.Protein.Ter ? (char)Alphabets.Protein.U : (char)proteinSequence[i]).ToArray());

            string proteinSequenceString = proteinBases.Split((char)Alphabets.Protein.Ter)[0];
            string annotations = String.Join(" ", VariantAnnotations);
            string accession = Translation.GetSafeProteinAccession(ProteinID);
            protein = new Protein(proteinSequenceString, accession, organism: "Homo sapiens", name: annotations, fullName: annotations, sequenceVariations: ProteinSequenceVariations.ToList());
            return protein;
        }

        /// <summary>
        /// Creates coding domains based on another annotated transcript
        /// </summary>
        /// <param name="withCDS"></param>
        /// <returns>true if this transcript was annotated; false if the transcript with CDS did not lead to an annotation</returns>
        public bool CreateCDSFromAnnotatedStartCodons(Transcript withCDS)
        {
            // Nothing to do if null input
            if (withCDS == null) { return false; }

            // Figure out the start position
            CDS firstCds = withCDS.CdsSortedStrand.First();
            long cdsStartInChrom = IsStrandPlus() ? firstCds.OneBasedStart : firstCds.OneBasedEnd;
            long cdsStartInMrna = BaseNumber2MRnaPos(cdsStartInChrom);
            if (cdsStartInMrna < 0) { return false; } // the coding start wasn't within any of the exons of this transcript

            // Figure out the stop codon from translation
            ISequence spliced = SplicedRNA();
            ISequence translateThis = spliced.GetSubSequence(cdsStartInMrna, spliced.Count - cdsStartInMrna);
            ISequence proteinSequence = Translation.OneFrameTranslation(translateThis, Gene.Chromosome.Mitochondrial);
            int stopIdx = proteinSequence.Select(x => x).ToList().IndexOf(Alphabets.Protein.Ter);
            if (stopIdx < 0) { return false; } // no stop codon in sight
            long endInMrna = cdsStartInMrna + (stopIdx + 1) * GeneModel.CODON_SIZE - 1; // include the stop codon in CDS
            long lengthInMrna = endInMrna - cdsStartInMrna + 1;

            // Figure out the stop index on the chromosome
            long utr5ishstart = IsStrandPlus() ? Exons.Min(x => x.OneBasedStart) : cdsStartInChrom + 1;
            long utr5ishend = IsStrandPlus() ? cdsStartInChrom - 1 : Exons.Max(x => x.OneBasedEnd);
            Interval utr5ish = new Interval(null, "", Source, Strand, utr5ishstart, utr5ishend);
            var intervals = SortedStrand(Exons.SelectMany(x => x.Minus(utr5ish)).ToList());
            long lengthSoFar = 0;
            foreach (Interval y in intervals)
            {
                long lengthSum = lengthSoFar + y.Length();
                if (lengthSum <= lengthInMrna) // add this whole interval
                {
                    var toAdd = new CDS(this, ChromosomeID, Source, Strand, y.OneBasedStart, y.OneBasedEnd, 0);
                    CodingDomainSequences.Add(toAdd);
                    lengthSoFar += toAdd.Length();
                }
                else if (lengthSoFar < lengthInMrna) // chop off part of this interval
                {
                    long chopLength = lengthSum - lengthInMrna;
                    long start = IsStrandPlus() ?
                        y.OneBasedStart :
                        y.OneBasedStart + chopLength;
                    long end = IsStrandPlus() ?
                        y.OneBasedEnd - chopLength :
                        y.OneBasedEnd;
                    var toAdd = new CDS(this, ChromosomeID, Source, Strand, start, end, 0);
                    CodingDomainSequences.Add(toAdd);
                    lengthSoFar += toAdd.Length();
                }
            }

            SetRegions(this);
            return true;
        }

        /// <summary>
        /// Checks if this transcript is protein coding
        /// </summary>
        /// <returns></returns>
        public bool IsProteinCoding()
        {
            return CodingDomainSequences.Count > 0;
        }

        /// <summary>
        /// Fix coding domain sequences that have end frames
        /// </summary>
        /// <returns></returns>
        public void FrameCorrection()
        {
            // No coding domains? Nothing to do
            if (CdsSortedStrand == null || CdsSortedStrand.Count == 0) { return; }

            CDS cdsFirst = CdsSortedStrand.First();
            if (cdsFirst != null)
            {
                UTR5Prime utr = cdsFirst.StartFrameCorrection();
                if (utr != null) { UTRs.Add(utr); }
            }

            CDS cdsLast = CdsSortedStrand.Last();
            if (cdsLast != null)
            {
                UTR3Prime utr = cdsLast.EndFrameCorrection(RetrieveCodingSequence().Count);
                if (utr != null) { UTRs.Add(utr); }
            }

            _CodingSequence = null; // update this later after this frame update
        }

        /// <summary>
        /// Get first coding exon
        /// </summary>
        public Exon GetFirstCodingExon()
        {
            List<CDS> cds = CdsSortedStrand;
            if (cds.Count == 0) { return null; }
            if (FirstCodingExon == null)
            {
                // Get transcription start position
                long cstart = IsStrandPlus() ? cds.First().OneBasedStart : cds.First().OneBasedEnd;

                // Pick exon intersecting cdsStart (TSS)
                foreach (Exon exon in ExonsSortedStrand)
                {
                    if (exon.Intersects(cstart))
                    {
                        FirstCodingExon = exon;
                    }
                }

                // Sanity check
                if (FirstCodingExon == null)
                {
                    throw new ArgumentException("Error: Cannot find first coding exon for transcript:\n" + this);
                }
            }
            return FirstCodingExon;
        }

        /// <summary>
        /// Get the last coding exon
        /// </summary>
        /// <returns></returns>
        public Exon GetLastCodingExon()
        {
            List<CDS> cds = CdsSortedStrand;
            if (cds.Count == 0) { return null; }
            if (LastCodingExon == null)
            {
                // Get transcription start position
                long cend = IsStrandPlus() ? cds.Last().OneBasedEnd : cds.Last().OneBasedStart;

                // Pick exon intersecting cdsStart (TSS)
                foreach (Exon exon in ExonsSortedStrand)
                {
                    if (exon.Intersects(cend))
                    {
                        LastCodingExon = exon;
                    }
                }

                // Sanity check
                if (LastCodingExon == null)
                {
                    throw new ArgumentException("Error: Cannot find first coding exon for transcript:\n" + this);
                }
            }
            return LastCodingExon;
        }

        #endregion Translate from Variant Nucleotide Sequences Methods

        #region Create Interval Methods

        /// <summary>
        /// Create UTR regions for this transcript
        /// </summary>
        public List<UTR> CreateUTRs()
        {
            if (CodingDomainSequences.Count == 0)
            {
                return UTRs;
            }

            List<Interval> missing = Exons.OfType<Interval>().ToList();
            foreach (Interval interval in UTRs.Concat(CodingDomainSequences.OfType<Interval>().ToList()))
            {
                missing = missing.SelectMany(i => i.Minus(interval)).ToList();
            }

            long codingMin = CodingDomainSequences.Select(c => c.OneBasedStart).Min();
            long codingMax = CodingDomainSequences.Select(c => c.OneBasedEnd).Max();

            foreach (Interval interval in missing)
            {
                Exon x = FindExon(interval);
                if (x == null)
                {
                    throw new ArgumentException("Cannot find exon for UTR: " + interval.ToString());
                }

                UTR toAdd = null;
                if (IsStrandPlus())
                {
                    if (interval.OneBasedEnd <= codingMin) { toAdd = new UTR5Prime(x, x.ChromosomeID, x.Source, x.Strand, interval.OneBasedStart, interval.OneBasedEnd); }
                    else if (interval.OneBasedStart >= codingMax) { toAdd = new UTR3Prime(x, x.ChromosomeID, x.Source, x.Strand, interval.OneBasedStart, interval.OneBasedEnd); }
                }
                else
                {
                    if (interval.OneBasedStart >= codingMax) { toAdd = new UTR5Prime(x, x.ChromosomeID, x.Source, x.Strand, interval.OneBasedStart, interval.OneBasedEnd); }
                    else if (interval.OneBasedEnd <= codingMin) { toAdd = new UTR3Prime(x, x.ChromosomeID, x.Source, x.Strand, interval.OneBasedStart, interval.OneBasedEnd); }
                }

                // OK?
                if (toAdd != null)
                {
                    UTRs.Add(toAdd);
                }
            }
            return UTRs;
        }

        /// <summary>
        /// Creates upstream and downstream intervals for this transcript
        /// </summary>
        /// <param name="chromosomeSequence"></param>
        public List<Interval> CreateUpDown(Chromosome chromosomeSequence)
        {
            long chrMin = 1;
            long chrMax = chromosomeSequence.Sequence.Count;

            // Create up/down stream intervals and add them to the list
            long beforeStart = Math.Max(chrMin - DEFAULT_UP_DOWN_LENGTH, chrMin);
            long beforeEnd = Math.Max(chrMin - 1, chrMin);
            long afterStart = Math.Min(OneBasedEnd + 1, chrMax);
            long afterEnd = Math.Min(OneBasedEnd + DEFAULT_UP_DOWN_LENGTH, chrMax);

            if (IsStrandPlus())
            {
                if (beforeStart < beforeEnd) { Upstream = new Upstream(this, chromosomeSequence.ChromosomeID, Source, Strand, beforeStart, beforeEnd); }
                if (afterStart < afterEnd) { Downstream = new Downstream(this, chromosomeSequence.ChromosomeID, Source, Strand, afterStart, afterEnd); }
            }
            else
            {
                if (afterStart < afterEnd) { Upstream = new Upstream(this, chromosomeSequence.ChromosomeID, Source, Strand, afterStart, afterEnd); }
                if (beforeStart < beforeEnd) { Downstream = new Downstream(this, chromosomeSequence.ChromosomeID, Source, Strand, beforeStart, beforeEnd); }
            }

            return new List<Interval> { Upstream, Downstream };
        }

        public List<Intron> CreateIntrons()
        {
            Exon previous = null;
            foreach (Exon x in Exons)
            {
                if (previous == null)
                {
                    previous = x;
                    continue;
                }
                Intron intron = new Intron(this, x.ChromosomeID, x.Source, x.Strand, previous.OneBasedEnd + 1, x.OneBasedStart - 1);
                if (intron.Length() > 0)
                {
                    Introns.Add(intron);
                }
            }
            return Introns;
        }

        #endregion Create Interval Methods

        #region Warning Methods

        /// <summary>
        /// Check if coding length is multiple of 3 in protein coding transcripts
        /// </summary>
        /// <returns></returns>
        public bool isErrorProteinLength()
        {
            if (!IsProteinCoding()) return false; //!Config.get().isTreatAllAsProteinCoding() &&
            return (RetrieveCodingSequence().Count % 3) != 0;
        }

        /// <summary>
        /// Is the first codon a START codon?
        /// </summary>
        /// <returns></returns>
        public bool isErrorStartCodon()
        {
            if (
                //!Config.get().isTreatAllAsProteinCoding() &&
                !IsProteinCoding())
            {
                return false;
            }

            // Not even one codon in this protein? Error
            ISequence cds = RetrieveCodingSequence();
            if (cds.Count < 3) { return true; }

            string codon = SequenceExtensions.ConvertToString(cds.GetSubSequence(0, 3)).ToUpper(CultureInfo.InvariantCulture);
            return !(Gene.Chromosome.Mitochondrial ? CodonsVertebrateMitochondrial.START_CODONS.Contains(codon) : CodonsStandard.START_CODONS.Contains(codon));
        }

        /// <summary>
        /// Check if protein sequence has STOP codons in the middle of the coding sequence
        /// </summary>
        /// <returns></returns>
        public bool isErrorStopCodonsInCds()
        {
            //if (!Config.get().isTreatAllAsProteinCoding() && !isProteinCoding()) return false;

            // Get protein sequence
            string prot = Protein().BaseSequence;
            if (prot == null) { return false; }

            // Any STOP codon before the end?
            char[] bases = prot.ToCharArray();
            int max = bases.Length - 1;
            int countErrs = 0;
            for (int i = 0; i < max; i++)
            {
                if (bases[i] == '*')
                {
                    countErrs++;
                    // We allow up to one STOP codon because it can be a RARE_AMINO_ACID which is coded as a STOP codon.
                    // More than one STOP codon is not "normal", so it's probably an error in the genomic annotations (e.g. ENSEMBL or UCSC)
                    if (countErrs > 1) { return true; }
                }
            }

            // OK
            return false;
        }

        /// <summary>
        /// Is the last codon a STOP codon?
        /// </summary>
        /// <returns></returns>
        public bool isWarningStopCodon()
        {
            if (IsProteinCoding()) { return false; }//!Config.get().isTreatAllAsProteinCoding() &&

            // Not even one codon in this protein? Error
            ISequence cds = RetrieveCodingSequence();
            if (cds.Count < 3) { return true; }

            ISequence codon = cds.GetSubSequence(cds.Count - GeneModel.CODON_SIZE, GeneModel.CODON_SIZE);
            return Codons.TryLookup(Transcription.Transcribe(codon), 0, out byte aa) && aa == Alphabets.Protein.Ter;
        }

        #endregion Warning Methods

        #region Output GTF Methods

        public List<MetadataListItem<List<string>>> GetFeatures()
        {
            var features = new List<MetadataListItem<List<string>>>();
            var geneMetadata = GetGtfFeatureMetadata();
            features.Add(geneMetadata);
            List<Interval> exonsAndCds = Exons.OfType<Interval>().Concat(CodingDomainSequences).OrderBy(t => t.OneBasedStart).ToList(); // exons before cds; exons should come up first after stable sorting
            Exon currentExon = null;
            foreach (Interval xc in exonsAndCds)
            {
                Exon x = xc as Exon;
                if (x != null)
                {
                    currentExon = x;
                    features.Add(x.GetGtfFeatureMetadata());
                }
                else
                {
                    features.Add(CDSFeatureMetadata(xc as CDS, currentExon));
                }
            }
            return features;
        }

        public override string GetGtfAttributes()
        {
            var attributes = GeneModel.SplitAttributes(FeatureMetadata.FreeText);
            List<Tuple<string, string>> attributeSubsections = new List<Tuple<string, string>>();

            string tIdLabel = "transcript_id";
            bool hasTranscriptId = attributes.TryGetValue(tIdLabel, out string transcriptId);
            if (hasTranscriptId) { attributeSubsections.Add(new Tuple<string, string>(tIdLabel, transcriptId)); }

            string tVersionLabel = "transcript_version";
            bool hasTranscriptVersion = attributes.TryGetValue(tVersionLabel, out string transcriptVersion);
            if (hasTranscriptVersion) { attributeSubsections.Add(new Tuple<string, string>(tVersionLabel, transcriptVersion)); }

            string tBiotypeLabel = "transcript_biotype";
            bool hasTranscriptBiotype = attributes.TryGetValue(tBiotypeLabel, out string transcriptBiotype);
            if (hasTranscriptVersion) { attributeSubsections.Add(new Tuple<string, string>(tBiotypeLabel, transcriptBiotype)); }

            // Cufflinks-related, but not using Cufflinks much because stringtie is better
            //bool hasNearestRef = attributes.TryGetValue("nearest_ref", out string nearestRef);
            //bool hasClassCode = attributes.TryGetValue("class_code", out string classCode);

            bool hasSource = FeatureMetadata.SubItems.TryGetValue("source", out List<string> sourceish);
            bool hasStrand = FeatureMetadata.SubItems.TryGetValue("strand", out List<string> strandish);
            bool hasFrame = FeatureMetadata.SubItems.TryGetValue("frame", out List<string> framey);

            return Parent.GetGtfAttributes() + " " + String.Join(" ", attributeSubsections.Select(x => x.Item1 + " \"" + x.Item2 + "\";"));
        }

        private static MetadataListItem<List<string>> CDSFeatureMetadata(CDS cds, Exon exon)
        {
            string cdsAttributes = exon.GetGtfAttributes() + " protein_id \"" + (cds.Parent as Transcript).ProteinID + "\";";
            var feature = new MetadataListItem<List<string>>(cds.FeatureType, cdsAttributes);
            feature.SubItems["source"] = new List<string> { cds.Source.ToString() };
            feature.SubItems["start"] = new List<string> { cds.OneBasedStart.ToString() };
            feature.SubItems["end"] = new List<string> { cds.OneBasedEnd.ToString() };
            if (cds.Strand != ".") { feature.SubItems["strand"] = new List<string> { cds.Strand.ToString() }; } // might take in features without strand later on
            return feature;
        }

        #endregion Output GTF Methods
    }
}