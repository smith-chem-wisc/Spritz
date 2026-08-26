using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace SpritzBackend
{
    /// <summary>
    /// The four values a run needs from a reference: which Ensembl release, which species directory, the
    /// common name, and the genome assembly.
    /// </summary>
    public sealed record SpritzReference(string Release, string Species, string Organism, string Genome);

    /// <summary>
    /// Turns the -r argument into a <see cref="SpritzReference"/>.
    ///
    /// Two forms are accepted. "release-116,homo_sapiens" names the release and the species and the
    /// assembly is looked up in genomes.csv. A whole genomes.csv line,
    /// "release-116,homo_sapiens,human,GRCh38", is also accepted so that saved commands keep working.
    /// </summary>
    public static class ReferenceString
    {
        public const string GenomesFileName = "genomes.csv";

        /// <summary>genomes.csv as installed beside the Spritz assembly.</summary>
        public static string InstalledGenomesPath() =>
            Path.Combine(AppContext.BaseDirectory, GenomesFileName);

        /// <summary>
        /// Resolves a reference. <paramref name="readGenomes"/> is only called for the two-field form, so
        /// the four-field form needs no genomes.csv present.
        /// </summary>
        public static SpritzReference Resolve(string reference, Func<IEnumerable<string>> readGenomes)
        {
            if (string.IsNullOrWhiteSpace(reference))
            {
                throw new SpritzException(
                    "Error: no reference specified. Pass a release and a species, e.g. " +
                    "\"release-116,homo_sapiens\", or a whole line from genomes.csv.");
            }

            string[] fields = reference.Split(',');
            for (int i = 0; i < fields.Length; i++)
            {
                if (string.IsNullOrWhiteSpace(fields[i]))
                {
                    throw new SpritzException($"Error: element {i + 1} of the reference string " +
                        $"\"{reference}\" is empty. Copy a whole line from genomes.csv, or pass just the " +
                        $"release and species, e.g. \"release-116,homo_sapiens\".");
                }
            }

            return fields.Length switch
            {
                2 => LookUp(ReleaseNumber(fields[0], reference), fields[1].Trim(), reference, readGenomes),
                4 => new SpritzReference(ReleaseNumber(fields[0], reference), fields[1].Trim(),
                        fields[2].Trim(), fields[3].Trim()),
                _ => throw new SpritzException($"Error: the reference string \"{reference}\" has " +
                        $"{fields.Length} comma-separated element(s). Pass two, a release and a species " +
                        $"such as \"release-116,homo_sapiens\", or a whole line from genomes.csv such as " +
                        $"\"release-116,homo_sapiens,human,GRCh38\"."),
            };
        }

        /// <summary>
        /// The release number, from either "release-116" or "116".
        /// </summary>
        private static string ReleaseNumber(string releaseField, string wholeReference)
        {
            const string prefix = "release-";
            string trimmed = releaseField.Trim();
            string number = trimmed.StartsWith(prefix, StringComparison.OrdinalIgnoreCase)
                ? trimmed[prefix.Length..]
                : trimmed;

            if (number.Length == 0 || !number.All(char.IsDigit))
            {
                throw new SpritzException($"Error: \"{releaseField}\" in the reference string " +
                    $"\"{wholeReference}\" is not an Ensembl release. Expected something like " +
                    $"\"release-116\" or \"116\".");
            }

            return number;
        }

        /// <summary>
        /// Finds the genome assembly and common name for a release and species in genomes.csv. The common
        /// name is carried through for the record; nothing in the workflow reads it.
        /// </summary>
        private static SpritzReference LookUp(string release, string species, string wholeReference,
            Func<IEnumerable<string>> readGenomes)
        {
            if (readGenomes is null)
            {
                throw new SpritzException($"Error: the reference string \"{wholeReference}\" needs " +
                    $"{GenomesFileName} to look up the genome assembly, and no source for it was given.");
            }

            int rows = 0;
            var releasesForSpecies = new List<string>();

            foreach (string line in readGenomes())
            {
                if (line.Length == 0 || line.StartsWith("#", StringComparison.Ordinal))
                {
                    continue;
                }

                string[] fields = line.Split(',');
                if (fields.Length != 4)
                {
                    continue;
                }

                rows++;
                bool sameSpecies = string.Equals(fields[1].Trim(), species, StringComparison.OrdinalIgnoreCase);
                if (!sameSpecies)
                {
                    continue;
                }

                string rowRelease = fields[0].Trim();
                if (rowRelease.StartsWith("release-", StringComparison.OrdinalIgnoreCase))
                {
                    rowRelease = rowRelease["release-".Length..];
                }

                if (rowRelease == release)
                {
                    return new SpritzReference(release, fields[1].Trim(), fields[2].Trim(), fields[3].Trim());
                }

                releasesForSpecies.Add(rowRelease);
            }

            if (rows == 0)
            {
                throw new SpritzException($"Error: {GenomesFileName} lists no references, so " +
                    $"\"{wholeReference}\" cannot be resolved. Regenerate it with the " +
                    $"--{SpritzOptionStrings.FetchGenomesLong} flag.");
            }

            if (releasesForSpecies.Count > 0)
            {
                throw new SpritzException($"Error: {GenomesFileName} has no release {release} for " +
                    $"\"{species}\". It has release(s) " +
                    $"{string.Join(", ", releasesForSpecies.Distinct().OrderBy(r => r.PadLeft(4, '0')))} " +
                    $"for that species.");
            }

            throw new SpritzException($"Error: {GenomesFileName} does not list the species \"{species}\" " +
                $"(searched {rows} references). Species are named as Ensembl names them, e.g. " +
                $"\"homo_sapiens\". Refresh the list with --{SpritzOptionStrings.FetchGenomesLong} if the " +
                $"species is newer than this copy.");
        }
    }
}
