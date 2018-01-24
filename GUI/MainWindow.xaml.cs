using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace GUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {

        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {

        }

        private void MenuItem_Help_Click(object sender, RoutedEventArgs e)
        {

        }

        private void MenuItem_Contact_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnAddFastq2Proteins_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new Fastq2ProteinsWF();
            if (dialog.ShowDialog()==true)
            {

            }
        }

        private void BtnAddLncRNADiscover_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new LncRNADiscoverWF();
            if (dialog.ShowDialog() == true)
            {

            }
        }

        private void RunWorkflowButton_Click(object sender, RoutedEventArgs e)
        {

        }

        private void DataGridCell_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {

        }

        private void BtnAddFASTA_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnClearFASTA_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnAddFastq_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnClearFastq_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnAddGeneSet_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnClearGeneSet_Click(object sender, RoutedEventArgs e)
        {

        }

        private void DataGridRow_Selected(object sender, RoutedEventArgs e)
        {

        }

        private void DataGridRow_Unselected(object sender, RoutedEventArgs e)
        {

        }
    }
}
