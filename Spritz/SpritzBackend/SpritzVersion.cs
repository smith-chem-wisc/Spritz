using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http;

namespace SpritzBackend
{
    public class SpritzVersion
    {
        public string NewestKnownVersion { get; private set; }
        public string NewestKnownVersionWithMsi { get; private set; }
        public bool IsMsiAvailableForUpdate { get; set; }

        public void GetVersionNumbersFromWeb()
        {
            // Attempt to get current MetaMorpheus version
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/Spritz/releases/latest").Result)
                {
                    var json = response.Content.ReadAsStringAsync().Result;
                    JObject deserialized = JObject.Parse(json);

                    // The API does not always answer with a release. A rate-limit or error reply is
                    // still valid JSON - {"message": "API rate limit exceeded for ...", ...} - and has
                    // no tag_name, so dereferencing it threw NullReferenceException. That surfaced as
                    // "Object reference not set to an instance of an object" both in CI and in the
                    // GUI's update dialog, which says nothing about what actually went wrong.
                    // Leaving the properties null lets the caller distinguish "no answer" from an
                    // answer it dislikes.
                    JToken tagName = deserialized["tag_name"];
                    if (tagName == null)
                    {
                        return;
                    }

                    NewestKnownVersion = tagName.ToString();
                    var assets = deserialized["assets"]?
                        .Select(b => b["name"]?.ToString())
                        .Where(name => name != null)
                        .ToList() ?? new List<string>();
                    bool containsMsi = assets.Contains("Spritz.msi");
                    if (!IsVersionLower(NewestKnownVersion))
                        IsMsiAvailableForUpdate = containsMsi;
                    if (containsMsi)
                        NewestKnownVersionWithMsi = NewestKnownVersion;
                }
            }
        }

        public static bool IsVersionLower(string checkVersionString)
        {
            var currV = GetVersionNumber(RunnerEngine.CurrentVersion);
            var checkVersion = GetVersionNumber(checkVersionString);
            return checkVersion.Item1 < currV.Item1 ||
                    (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 < currV.Item2) ||
                    (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 == currV.Item2 && checkVersion.Item3 <= currV.Item3);
        }

        public static (int, int, int) GetVersionNumber(string VersionNode)
        {
            try
            {
                var split = VersionNode.Split('.');

                return (int.Parse(split[0]), int.Parse(split[1]), int.Parse(split[2]));
            }
            catch (FormatException)
            {
                return (0, 0, 0);
            }
        }
    }
}