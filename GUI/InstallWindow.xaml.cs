using CMD;
using System.Diagnostics;
using System.Windows;
using System.Windows.Navigation;
using ToolWrapperLayer;

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
                TxbkInstall.Text = "The Windows Subsystem for Linux has not been enabled. Please see link below for more details.";
            }
            else
            {
                TxbkInstall.Text = "Time to install the required packages in the Windows Subsystem for Linux!";
                BtnInstall.IsEnabled = true;
                BtnAlreadyInstalled.IsEnabled = true;
            }
        }

        private void BtnInstall_Click(object sender, RoutedEventArgs e)
        {
            Spritz.Main(new[] { "CMD.exe", "-c", "setup" });
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