using Newtonsoft.Json.Linq;
using System;
using System.Linq;
using System.Net.Http;

namespace SpritzBackend
{
    public class SpritzVersion
    {
        public string NewestKnownVersion { get; private set; }
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
                    NewestKnownVersion = deserialized["tag_name"].ToString();
                    if (!IsVersionLower(RunnerEngine.CurrentVersion, NewestKnownVersion))
                    {
                        var assets = deserialized["assets"].Select(b => b["name"].ToString()).ToList();
                        IsMsiAvailableForUpdate = assets.Contains("Spritz.msi");
                    }
                }
            }
        }

        public static bool IsVersionLower(string currVersion, string checkVersionString)
        {
            var currV = GetVersionNumber(currVersion);
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