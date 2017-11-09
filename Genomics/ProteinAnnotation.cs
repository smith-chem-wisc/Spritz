using Bio;
using Bio.Extensions;
using Bio.Algorithms.Translation;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Proteogenomics
{
    public static class ProteinAnnotation
    {
        public static void Annotate(TranscriptPossiblyWithVariants transcript)
        {
            Annotate(new List<TranscriptPossiblyWithVariants> { transcript });
        }

        // get the SAV notation X#X
        // or for indels use X#XXX or XXX#X where # is the start
        // get the SNV location 1:100000
        // get the codon change Gcc/Acc
        public static void Annotate(List<TranscriptPossiblyWithVariants> transcripts)
        {
            foreach (TranscriptPossiblyWithVariants t in transcripts)
            {
                ISequence referenceTranscriptSequence = new Sequence(Alphabets.DNA, String.Join("", t.GetExonsUsedInDerivation().Select(x => SequenceExtensions.ConvertToString(x.Sequence))));
                ISequence variantTranscriptSequence = new Sequence(Alphabets.DNA, t.Sequence);
                int[] indices = t.GetExonsUsedInDerivation().SelectMany(x => Enumerable.Range((int)x.OneBasedStart, (int)(x.OneBasedEnd - x.OneBasedStart + 1))).ToArray();
                referenceTranscriptSequence = t.Transcript.Strand == "+" ? referenceTranscriptSequence : referenceTranscriptSequence.GetReverseComplementedSequence();
                variantTranscriptSequence = t.Transcript.Strand == "+" ? variantTranscriptSequence : variantTranscriptSequence.GetReverseComplementedSequence();
                indices = t.Transcript.Strand == "+" ? indices : indices.Reverse().ToArray();

                for (long i = t.ZeroBasedCodingStart; i + 2 < referenceTranscriptSequence.Count; i += 3)
                {
                    long oneBasedTranscriptStart = t.Transcript.Exons.Min(x => x.OneBasedStart);
                    long[] codonIndices = new long[] { indices[i], indices[i + 1], indices[i + 2] };
                    Variant v = t.Variants.FirstOrDefault(vv => codonIndices.Contains(vv.OneBasedStart));
                    if (v != null)
                    {
                        long j = i - t.ZeroBasedCodingStart;
                        long aminoAcidPosition = i / 3 + 1;
                        Codons.TryLookup(Transcription.GetRnaComplement(referenceTranscriptSequence[i]), 
                            Transcription.GetRnaComplement(referenceTranscriptSequence[i + 1]), 
                            Transcription.GetRnaComplement(referenceTranscriptSequence[i + 2]), 
                            out byte originalAminoAcid);
                        Codons.TryLookup(Transcription.GetRnaComplement(variantTranscriptSequence[j]),
                            Transcription.GetRnaComplement(variantTranscriptSequence[j + 1]),
                            Transcription.GetRnaComplement(variantTranscriptSequence[j + 2]),
                            out byte newAminoAcid);
                        if (v.ReferenceAllele.Length == v.AlternateAllele.Length) // single amino acid variation
                        {
                            v.Synonymous = originalAminoAcid == newAminoAcid;
                            v.Annotation = "pep:sav " +
                                char.ToUpperInvariant((char)originalAminoAcid) + aminoAcidPosition.ToString() + char.ToUpperInvariant((char)newAminoAcid) + " " +
                                v.Chr + ":" + v.OneBasedStart + " " +
                                String.Join("",
                                    char.ToUpperInvariant((char)referenceTranscriptSequence[i]),
                                    char.ToUpperInvariant((char)referenceTranscriptSequence[i + 1]),
                                    char.ToUpperInvariant((char)referenceTranscriptSequence[i + 2]))
                                + "/" +
                                String.Join("",
                                    char.ToUpperInvariant((char)variantTranscriptSequence[j]),
                                    char.ToUpperInvariant((char)variantTranscriptSequence[j + 1]),
                                    char.ToUpperInvariant((char)variantTranscriptSequence[j + 2]));
                        }

                        // TODO: adjust indices to keep stepping across the right portion of the variant transcript
                        else
                        {
                            v.Synonymous = originalAminoAcid == newAminoAcid && char.ToUpperInvariant((char)originalAminoAcid) == '*';
                            bool insertion = v.ReferenceAllele.Length < v.AlternateAllele.Length;
                            bool deletion = !insertion;
                            string longerAllele = insertion ? v.AlternateAllele : v.ReferenceAllele;
                            string shorterAllele = insertion ? v.ReferenceAllele : v.AlternateAllele;
                            bool frameshift = (longerAllele.Length - shorterAllele.Length) % 3 != 0;
                            if (frameshift)
                            {
                                v.Annotation = "pep:" + (insertion ? "frameShiftInsertion " : "frameShiftDeletion ") +
                                    originalAminoAcid + aminoAcidPosition.ToString() +
                                    (insertion ? "frameShiftInsertion " : "frameShiftDeletion ") +
                                    v.Chr + ":" + v.OneBasedStart;
                            }
                            else
                            { 
                                string originalSequence = new string(new char[] { char.ToUpperInvariant((char)originalAminoAcid) });
                                string variantSequence = new string(new char[] { char.ToUpperInvariant((char)newAminoAcid) });
                                for (long ii = i + 3; ii + 2 < referenceTranscriptSequence.Count && (ii - i) / 3 < ((longerAllele.Length - shorterAllele.Length) / 3) + 1; ii += 3)
                                {
                                    long jj = ii - t.ZeroBasedCodingStart;
                                    if (insertion && Codons.TryLookup(
                                            Transcription.GetRnaComplement(variantTranscriptSequence[jj]), 
                                            Transcription.GetRnaComplement(variantTranscriptSequence[jj + 1]), 
                                            Transcription.GetRnaComplement(variantTranscriptSequence[jj + 2]), 
                                            out byte newAminoAcid2))
                                        variantSequence += newAminoAcid2;
                                    if (deletion && Codons.TryLookup(Transcription.GetRnaComplement(referenceTranscriptSequence[ii]), 
                                            Transcription.GetRnaComplement(referenceTranscriptSequence[ii + 1]), 
                                            Transcription.GetRnaComplement(referenceTranscriptSequence[ii + 2]), 
                                            out byte originalAminoAcid2))
                                        originalSequence += originalAminoAcid2;
                                }
                                v.Annotation = "pep:" + (insertion ? "inFrameInsertion " : "inFrameDeletion ") + 
                                    originalSequence + aminoAcidPosition.ToString() + variantSequence + " " +
                                    v.Chr + ":" + v.OneBasedStart;
                            }
                        }
                    }
                }
                t.ProteinAnnotation = String.Join(" ", t.Variants.Select(v => v.Annotation));
            }
        }

        /// <summary>
        /// Transfers likely modifications from a list of proteins to another based on sequence similarity. Returns a list of new objects.
        /// </summary>
        /// <param name="destination"></param>
        /// <param name="source"></param>
        /// <returns></returns>
        public static List<Protein> TransferModifications(List<Protein> source, List<Protein> destination)
        {
            List<Protein> newProteins = new List<Protein>();
            Dictionary<string, Protein> dictDestination = destination.ToDictionary(p => p.BaseSequence, p => p);
            Dictionary<string, Protein> dictSource = source.ToDictionary(p => p.BaseSequence, p => p);

            List<string> commonSeqs = dictDestination.Keys.Intersect(dictDestination.Keys).ToList();
            foreach (string seq in commonSeqs)
            {
                Protein source_protein = dictSource[seq];
                Protein destination_protein = dictDestination[seq];
                newProteins.Add(
                    new Protein(
                        seq,
                        destination_protein.Accession,
                        gene_names: source_protein.GeneNames.ToList(),
                        oneBasedModifications: source_protein.OneBasedPossibleLocalizedModifications,
                        proteolysisProducts: source_protein.ProteolysisProducts.ToList(),
                        name: destination_protein.Name,
                        full_name: destination_protein.FullName,
                        isDecoy: destination_protein.IsDecoy,
                        isContaminant: destination_protein.IsContaminant,
                        databaseReferences: source_protein.DatabaseReferences.ToList(),
                        sequenceVariations: source_protein.SequenceVariations.ToList()));
            }

            List<string> destinationOnly = dictDestination.Keys.Except(dictSource.Keys).ToList();
            List<string> sourceOnly = dictSource.Keys.Except(dictDestination.Keys).ToList();
            return newProteins;
        }
    }
}
