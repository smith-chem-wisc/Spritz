using Newtonsoft.Json.Linq;
using SpritzBackend;
using System;
using System.Diagnostics;
using System.IO;
using System.Net.Http;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

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
            lbl.Text = "A newer version: " + Version.NewestKnownVersion + " is available!";
            ReleaseHandler();
        }

        public static SpritzVersion Version { get; private set; } = new();

        private void InstallerClicked(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
            HttpClient client = new();

            var uri = new Uri(@"https://github.com/smith-chem-wisc/Spritz/releases/download/" + Version.NewestKnownVersion + @"/Spritz.msi");

            Exception exception = null;
            try
            {
                // download and start the installer
                var tempDownloadLocation = Path.Combine(Path.GetTempPath(), "Spritz.msi");
                var urlResponse = Task.Run(() => client.GetAsync(uri)).Result;
                using (FileStream stream = new(tempDownloadLocation, FileMode.CreateNew))
                {
                    Task.Run(() => urlResponse.Content.CopyToAsync(stream)).Wait();
                }

                // start the installer
                Process p = new();
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

        private void ReleaseHandler()
        {
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/Spritz/releases").Result)
                {
                    string json = response.Content.ReadAsStringAsync().Result;
                    JArray GitArray = JArray.Parse(json);
                    StringBuilder allVersionsText = new();
                    foreach (JObject obj in GitArray.Children<JObject>())
                    {
                        if (SpritzVersion.IsVersionLower(obj.SelectToken("tag_name").ToString()))
                            break;
                        string body = new MarkdownSharp.Markdown().Transform(obj.SelectToken("body").ToString());
                        allVersionsText.AppendLine("<font face=\"Arial\" size=2>");
                        allVersionsText.AppendLine("<h3>" + obj.SelectToken("tag_name").ToString() + "</h3>");
                        allVersionsText.AppendLine(body);
                        allVersionsText.AppendLine();
                        if (!Version.IsMsiAvailableForUpdate)
                        {
                            allVersionsText.AppendLine("Warning: A new version of Spritz was detected, but the installer file hasn't been uploaded yet. Please try again in a few minutes.");
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
            Version.GetVersionNumbersFromWeb();
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
    }
}