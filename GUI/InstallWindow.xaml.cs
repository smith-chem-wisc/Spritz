using System;
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
            Process proc = new Process();
            proc.StartInfo.FileName = "CMD.exe";
            proc.StartInfo.Arguments = "-c setup";
            proc.StartInfo.CreateNoWindow = true;
            proc.StartInfo.UseShellExecute = false; // don't fire up a shell for the CMD.exe process
            proc.Start();
            proc.WaitForExit();
            DialogResult = true;

            bool installationSuccessful = WrapperUtility.CheckBashSetup() && WrapperUtility.CheckToolSetup(Environment.CurrentDirectory);
            MessageBox.Show("Installation of tools " + (installationSuccessful ? "was" : "was not") + " successful.", "Spritz Installation", MessageBoxButton.OK);
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