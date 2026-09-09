using System;
using System.IO;
using System.Linq;

namespace SpritzModifications
{
    /// <summary>
    /// Rewrites the FASTA headers mzLib writes into the UniProt-style form Philosopher/FragPipe parses.
    ///
    /// mzLib writes ">mz|{accession}|{name} {fullName} OS={organism} [OX={taxid} ]GN={gene}", and for a
    /// variant entry both name and fullName carry a " variant:" tag whose payload is the original raw VCF
    /// line - pipes, commas and an unescaped '>' included. Philosopher's documented requirement is that the
    /// description "should not contain any commas or special characters", so the tag is dropped here and the
    /// per-variant detail is carried by the sibling .vardesc.tsv instead, keyed on the same accession.
    /// </summary>
    public static class FastaHeaders
    {
        public const string VariantTag = " variant:";
        private const string DatabasePrefix = "mz|";
        private const string PhilosopherDatabasePrefix = "sp|";
        private const string DecoyPrefix = "rev_mz|";
        private const string PhilosopherDecoyPrefix = "rev_sp|";

        /// <summary>
        /// Converts one FASTA header line (with or without its leading '>') to the Philosopher form.
        /// Lines that are not headers are returned unchanged, so this is safe to map over a whole file.
        /// Idempotent: a header that has already been converted is returned unchanged.
        /// </summary>
        public static string ToPhilosopherHeader(string line)
        {
            if (line == null || !line.StartsWith(">", StringComparison.Ordinal))
            {
                return line;
            }

            string header = line.Substring(1);
            header = StripVariantTag(header);

            // Ensembl identifiers arrive as "transcript:ENST...", which is not how anyone refers to them.
            header = header.Replace("transcript:", "");

            if (header.StartsWith(DecoyPrefix, StringComparison.Ordinal))
            {
                header = PhilosopherDecoyPrefix + header.Substring(DecoyPrefix.Length);
            }
            else if (header.StartsWith(DatabasePrefix, StringComparison.Ordinal))
            {
                header = PhilosopherDatabasePrefix + header.Substring(DatabasePrefix.Length);
            }

            return ">" + CollapseSpaces(header);
        }

        /// <summary>
        /// Removes the " variant:" payload, which runs from the first tag up to the OS= field that follows
        /// the name and fullName. Leaves a header without a tag untouched.
        /// </summary>
        private static string StripVariantTag(string header)
        {
            int tagStart = header.IndexOf(VariantTag, StringComparison.Ordinal);
            if (tagStart < 0)
            {
                return header;
            }

            // OS= is the first field after fullName, so it terminates the payload. Search from the tag so a
            // stray "OS=" inside an accession cannot be mistaken for the terminator.
            int organism = header.IndexOf(" OS=", tagStart, StringComparison.Ordinal);
            return organism < 0 ? header.Substring(0, tagStart) : header.Substring(0, tagStart) + header.Substring(organism);
        }

        /// <summary>
        /// Spritz XML entries usually have no fullName, which leaves a double space where it would have been.
        /// Applies to every header, not just ones that carried a variant tag.
        /// </summary>
        private static string CollapseSpaces(string text)
        {
            return string.Join(" ", text.Split(' ').Where(t => t.Length > 0));
        }

        /// <summary>
        /// Sibling path for the Philosopher-ready copy: foo.fasta -> foo.fragpipe.fasta.
        /// </summary>
        public static string PhilosopherFastaPath(string fastaPath)
        {
            return Path.Combine(
                Path.GetDirectoryName(fastaPath) ?? "",
                Path.GetFileNameWithoutExtension(fastaPath) + ".fragpipe.fasta");
        }

        /// <summary>
        /// Writes a Philosopher-ready copy of a FASTA file. Sequences are copied verbatim; only headers move.
        /// </summary>
        public static void WritePhilosopherFasta(string sourceFasta, string destinationFasta)
        {
            File.WriteAllLines(destinationFasta, File.ReadLines(sourceFasta).Select(ToPhilosopherHeader));
        }
    }
}
