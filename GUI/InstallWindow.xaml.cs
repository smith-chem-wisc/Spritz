using System.Windows;
using CMD;

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
        }

        private void BtnInstall_Click(object sender, RoutedEventArgs e)
        {
            Spritz.Main(new string[] { "CMD.exe", "-c", "setup" });
            DialogResult = true;
            //return;
        }

        private void BtnAlreadyInstalled_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}
