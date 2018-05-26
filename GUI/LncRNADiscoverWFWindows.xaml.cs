using System.Windows;
using WorkflowLayer;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for LncRNADiscoverWF.xaml
    /// </summary>
    public partial class LncRNADiscoverFlowWindows : Window
    {
        public LncRNADiscoverFlowWindows()
        {
            InitializeComponent();

            TheTask = new LncRNADiscoveryFlow();
        }

        internal LncRNADiscoveryFlow TheTask { get; private set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            //LncRNADiscoveryParameters parameters = new LncRNADiscoveryParameters();

            DialogResult = true;
        }
    }
}