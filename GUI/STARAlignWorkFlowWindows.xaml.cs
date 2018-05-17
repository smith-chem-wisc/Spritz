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
using System.Collections.ObjectModel;
using System.IO;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for STARAlignWorkFlowWindows.xaml
    /// </summary>
    public partial class STARAlignWorkFlowWindows : Window
    {
        public STARAlignWorkFlowWindows()
        {
            InitializeComponent();
            mainWindow = (MainWindow)Application.Current.MainWindow;
            TheTask = new STARAlignmentFlow();
            TheTask.Parameters = new STARAlignmentParameters();
            TheTask.SpritzParameters = TheTask.Parameters;
            UpdateFieldsFromTask(TheTask);
        }

        internal STARAlignmentFlow TheTask { get; private set; }

        private MainWindow mainWindow { get; set; }

        protected void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        protected void saveButton_Click(object sender, RoutedEventArgs e)
        {
            STARAlignmentParameters parametersToSave = new STARAlignmentParameters();
            parametersToSave.AnalysisDirectory = txtAnalysisDirectory.Text;
            parametersToSave.Threads = int.Parse(txtThreads.Text);
            parametersToSave.StrandSpecific = ckbStrandSpecific.IsChecked.Value;
            parametersToSave.OverWriteStarAlignment = ckbOverWriteStarAlignment.IsChecked.Value;
            parametersToSave.GenomeStarIndexDirectory = txtGenomeStarIndexDirectory.Text;
            parametersToSave.ReorderedFasta = txtReorderedFasta.Text;
            parametersToSave.EnsemblKnownSitesPath = txtEnsemblKnownSitesPath.Text;
            parametersToSave.UseReadSubset = ckbReadSubset.IsChecked.Value;
            parametersToSave.ReadSubset = int.Parse(txtReadSubset.Text);

            //Pass file Parameters from MainWindow. Create new function for this part.
            var genomeFastaDataGrids = (ObservableCollection<GenomeFastaDataGrid>)mainWindow.dataGridFASTA.DataContext;
            parametersToSave.ReorderedFasta = genomeFastaDataGrids.First().FilePath;
            var geneSetCollection = (ObservableCollection<GeneSetDataGrid>)mainWindow.dataGridGeneSet.DataContext;
            parametersToSave.GeneModelGtfOrGff = geneSetCollection.First().FilePath;
            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)mainWindow.dataGridRnaSeqFastq.DataContext;
            parametersToSave.Fastqs = new List<string[]> { rnaSeqFastqCollection.Select(p => p.FilePath).ToArray() };

            parametersToSave.AnalysisDirectory = System.IO.Path.Combine(mainWindow.OutputFolderTextBox.Text);
            TheTask.Parameters = parametersToSave;
            TheTask.SpritzParameters = TheTask.Parameters;
            DialogResult = true;
        }

        private void UpdateFieldsFromTask(STARAlignmentFlow sTARAlignmentFlow)
        {
            txtAnalysisDirectory.Text = sTARAlignmentFlow.Parameters.AnalysisDirectory;
            txtThreads.Text = sTARAlignmentFlow.Parameters.Threads.ToString();
            ckbStrandSpecific.IsChecked = sTARAlignmentFlow.Parameters.StrandSpecific;
            ckbInferStrandSpecificity.IsChecked = sTARAlignmentFlow.Parameters.InferStrandSpecificity;
            ckbOverWriteStarAlignment.IsChecked = sTARAlignmentFlow.Parameters.OverWriteStarAlignment;
            txtGenomeStarIndexDirectory.Text = sTARAlignmentFlow.Parameters.GenomeStarIndexDirectory;
            txtReorderedFasta.Text = sTARAlignmentFlow.Parameters.ReorderedFasta;
            txtEnsemblKnownSitesPath.Text = sTARAlignmentFlow.Parameters.EnsemblKnownSitesPath;
            ckbReadSubset.IsChecked = sTARAlignmentFlow.Parameters.UseReadSubset;
            txtReadSubset.Text = sTARAlignmentFlow.Parameters.ReadSubset.ToString();

        }
    }
}
