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
using System.Windows.Shapes;
using WorkflowLayer;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for Fastq2ProteinsWF.xaml
    /// </summary>
    public partial class Fastq2ProteinsWFWindows : Window
    {
        public Fastq2ProteinsWFWindows()
        {
            InitializeComponent();

            TheTask = new Fastq2ProteinsWF();
        }

        internal Fastq2ProteinsWF TheTask { get; private set; }


        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {

        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            RnaSeqAlignParameters rnaSeqAlignParameters = new RnaSeqAlignParameters();

            DialogResult = true;
        }
    }
}
