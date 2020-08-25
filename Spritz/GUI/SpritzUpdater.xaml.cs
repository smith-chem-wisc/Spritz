using Newtonsoft.Json.Linq;
using System;
using System.Diagnostics;
using System.Net;
using System.Net.Http;
using System.Text;
using System.Windows;
using System.Linq;
using System.IO;

namespace Spritz
{
    /// <summary>
    /// Interaction logic for SpritzUpdater.xaml
    /// </summary>
    public partial class SpritzUpdater : Window
    {
        public SpritzUpdater()
        {
            InitializeComponent();
            lbl.Text = "A newer version: " + NewestKnownVersion + " is available!";
            ReleaseHandler();
        }

        public static string NewestKnownVersion { get; private set; }
        private static bool IsMsiAvailableForUpdate { get; set; }

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

        private void InstallerClicked(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
            using (var client = new WebClient())
            {
                var uri = new Uri(@"https://github.com/smith-chem-wisc/Spritz/releases/download/" + NewestKnownVersion + @"/SpritzInstaller.msi");

                Exception exception = null;
                try
                {
                    // download and start the installer
                    var tempDownloadLocation = Path.Combine(System.IO.Path.GetTempPath(), "SpritzInstaller.msi");
                    client.DownloadFile(uri, tempDownloadLocation);
                    var p = new Process();
                    p.StartInfo = new ProcessStartInfo()
                    {
                        UseShellExecute = true,
                        FileName = tempDownloadLocation
                    };
                    p.Start();
                }
                catch (Exception ex)
                {
                    exception = ex;
                    MessageBox.Show(ex.Message);
                }

                if (exception == null)
                {
                    // close Spritz if the installer was started successfully
                    Application.Current.Shutdown();
                }
            }
        }

        private void ReleaseHandler()
        {
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/Spritz/releases").Result)
                {
                    string json = response.Content.ReadAsStringAsync().Result;
                    JArray GitArray = JArray.Parse(json);
                    var currV = GetVersionNumber(MainWindow.CurrentVersion);
                    StringBuilder allVersionsText = new StringBuilder();
                    foreach (JObject obj in GitArray.Children<JObject>())
                    {
                        if (IsVersionLower(obj.SelectToken("tag_name").ToString()))
                            break;
                        string body = new MarkdownSharp.Markdown().Transform(obj.SelectToken("body").ToString());
                        allVersionsText.AppendLine("<font face=\"Arial\" size=2>");
                        allVersionsText.AppendLine("<h3>" + obj.SelectToken("tag_name").ToString() + "</h3>");
                        allVersionsText.AppendLine(body);
                        allVersionsText.AppendLine();
                        if (!IsMsiAvailableForUpdate)
                        {
                            allVersionsText.AppendLine("Warning: A new version of MetaMorpheus was detected, but the installer file hasn't been uploaded yet. Please try again in a few minutes.");
                            Bt_DownloadAndRunInstaller.IsEnabled = false;
                        }
                        allVersionsText.AppendLine("</font>");
                    }
                    releases.NavigateToString(allVersionsText.ToString());
                    releases.Navigating += Releases_Navigating;
                }
            }
        }

        public static void GetVersionNumbersFromWeb()
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
                    if (!IsVersionLower(NewestKnownVersion))
                    {
                        var assets = deserialized["assets"].Select(b => b["name"].ToString()).ToList();
                        IsMsiAvailableForUpdate =  assets.Contains("SpritzInstaller.msi");
                    }
                }
            }
        }

        public void Releases_Navigating(object sender, System.Windows.Navigation.NavigatingCancelEventArgs e)
        {
            e.Cancel = true; //cancel the current event
            Process.Start(e.Uri.ToString()); //this opens the URL in the user's default browser
        }

        private void NoClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        public static bool IsVersionLower(string checkVersionString)
        {
            var currV = GetVersionNumber(MainWindow.CurrentVersion);
            var checkVersion = GetVersionNumber(checkVersionString);
            return checkVersion.Item1 < currV.Item1 ||
                    (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 < currV.Item2) ||
                    (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 == currV.Item2 && checkVersion.Item3 <= currV.Item3);
        }
    }
}