using System.Windows;
using CMD;
using ToolWrapperLayer;
using System.Diagnostics;
using System.Windows.Navigation;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for InstallWindow.xaml
    /// </summary>
    public partial class InstallWindow : Window
    {
        public InstallWindow()
        {
            InitializeComponent();
            CheckBashSetup();
        }

        private void CheckBashSetup()
        {
            if (!WrapperUtility.CheckBashSetup())
            {
                TxbkInstall.Text = "The Windows Subsystem for Windows has not been enabled. Please see open link below for more details.";
                
            }
            else
            {
                TxbkInstall.Text = "Please install all the required packages!";
                BtnInstall.IsEnabled = true;
                BtnAlreadyInstalled.IsEnabled = true;
            }
        }

        private void BtnInstall_Click(object sender, RoutedEventArgs e)
        {
            Spritz.Main(new string[] { "CMD.exe", "-c", "setup" });
            DialogResult = true;
        }

        private void BtnAlreadyInstalled_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Url_Click(object sender, RequestNavigateEventArgs e)
        {
            Process.Start(new ProcessStartInfo(e.Uri.AbsoluteUri));
            e.Handled = true;
        }
    }
}
