using System;
using System.Net.Http;

namespace SpritzBackend
{
    /// <summary>
    /// Whether Ensembl publishes known variant sites for a species, which decides whether GATK base
    /// recalibration can download them or has to bootstrap its own.
    ///
    /// Worth checking rather than assuming: at release 116 only 19 of the 359 species with a gene
    /// model have a variation/vcf/ directory. Saccharomyces cerevisiae had one at release 96 and
    /// does not now, so the answer moves between releases for the same species.
    /// </summary>
    public static class EnsemblVariation
    {
        public const string Auto = "auto";
        public const string Ensembl = "ensembl";
        public const string Bootstrap = "bootstrap";

        public static bool IsKnownMode(string mode) =>
            string.Equals(mode, Auto, StringComparison.OrdinalIgnoreCase) ||
            string.Equals(mode, Ensembl, StringComparison.OrdinalIgnoreCase) ||
            string.Equals(mode, Bootstrap, StringComparison.OrdinalIgnoreCase);

        /// <summary>
        /// Both filenames download_ensembl_vcf can fetch. The rule tries the capitalised one first
        /// and falls back to the lowercase one; only the lowercase form resolves for any species
        /// today, so it is checked first here. Order does not matter to the answer - this is a
        /// disjunction - but it is not the rule's order.
        /// </summary>
        public static string[] CandidateUrls(string release, string speciesLower)
        {
            string directory = $"https://ftp.ensembl.org/pub/release-{release}/variation/vcf/{speciesLower}";
            string capitalised = char.ToUpperInvariant(speciesLower[0]) + speciesLower[1..];
            return new[] { $"{directory}/{speciesLower}.vcf.gz", $"{directory}/{capitalised}.vcf.gz" };
        }

        /// <summary>
        /// A HEAD against each candidate. A network failure is reported as "published" rather than
        /// "absent", so an offline run keeps the behaviour it had before this check existed and
        /// fails in download_ensembl_vcf with its own message, instead of silently switching a run
        /// to a different recalibration strategy on the strength of a dropped connection.
        /// </summary>
        public static bool Published(string release, string speciesLower)
        {
            try
            {
                using var client = new HttpClient { Timeout = TimeSpan.FromSeconds(20) };
                foreach (string url in CandidateUrls(release, speciesLower))
                {
                    using var request = new HttpRequestMessage(HttpMethod.Head, url);
                    using HttpResponseMessage response = client.Send(request);
                    if (response.IsSuccessStatusCode)
                    {
                        return true;
                    }
                    // Only "not there" counts as absent. A 403 from a proxy, a 429, a 502 or a 405
                    // says nothing about whether Ensembl publishes variation, and treating it as
                    // absence would switch the run's recalibration strategy on the strength of an
                    // unrelated failure - the thing the exception handler below exists to avoid.
                    if (response.StatusCode != System.Net.HttpStatusCode.NotFound
                        && response.StatusCode != System.Net.HttpStatusCode.Gone)
                    {
                        return true;
                    }
                }
                return false;
            }
            catch (Exception)
            {
                return true;
            }
        }
    }
}
