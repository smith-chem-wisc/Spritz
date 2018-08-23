using CMD;
using System.Collections.ObjectModel;
using System.Linq;
using System.Windows;
using WorkflowLayer;
using System;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for TransferModificationsFlowWindows.xaml
    /// </summary>
    public partial class WorkFlowWindow : Window
    {
        public WorkFlowWindow()
        {
            InitializeComponent();
            PopulateChoices();
            Options = new Options();
            mainWindow = (MainWindow)Application.Current.MainWindow;
            UpdateFieldsFromTask();
        }

        private MainWindow mainWindow { get; set; }

        public Options Options { get; set; }

        protected void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        protected void saveButton_Click(object sender, RoutedEventArgs e)
        {
            int i = CbxWorkFlowType.SelectedIndex;
            switch (i)
            {
                case 0:
                    Options.Command = SampleSpecificProteinDBFlow.Command;
                    break;
                case 1:
                    Options.Command = LncRNADiscoveryFlow.Command;
                    break;
                case 2:
                    Options.Command = TranscriptQuantificationFlow.Command;
                    break;
                case 3:
                    //Options.Command = STARAlignmentFlow.Command;
                    break;
                case 4:
                    Options.Command = GeneFusionDiscoveryFlow.Command;
                    break;
                default:
                    break;
            }
            
            Options.Threads = int.Parse(txtThreads.Text);
            Options.AnalysisDirectory = txtAnalysisDirectory.Text;
            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)mainWindow.dataGridRnaSeqFastq.DataContext;
            Options.Fastq1 = rnaSeqFastqCollection.Select(p => p.FilePath).ToList()[0];
            Options.Fastq2 = rnaSeqFastqCollection.Select(p => p.FilePath).ToList()[1];

            Options.GenomeStarIndexDirectory = txtGenomeStarIndexDirectory.Text;
            Options.StrandSpecific = ckbStrandSpecific.IsChecked.Value;
            Options.UniProtXml = txtUniProtProteinXml.Text;
            Options.SraAccession = txtSraAccession.Text;
            Options.SpritzDirectory = txtSpritzDirecory.Text;
            Options.ReferenceVcf = txtDbsnpVcfReference.Text;
            Options.Reference = txtReference.Text;
            Options.ProteinFastaPath = txtProteinFasta.Text;
            Options.OverwriteStarAlignments = ckbOverWriteStarAlignment.IsChecked.Value;
            Options.GenomeStarIndexDirectory = txtGenomeStarIndexDirectory.Text;
            Options.InferStrandSpecificity = ckbInferStrandedness.IsChecked.Value;

            DialogResult = true;
        }

        private void UpdateFieldsFromTask()
        {
            CbxWorkFlowType.SelectedIndex = 0;
            txtAnalysisDirectory.Text = Options.AnalysisDirectory;
            txtThreads.Text = Options.Threads.ToString();
            ckbStrandSpecific.IsChecked = Options.StrandSpecific;
            ckbStrandSpecific.IsChecked = Options.InferStrandSpecificity;
            ckbOverWriteStarAlignment.IsChecked = Options.OverwriteStarAlignments;
            txtGenomeStarIndexDirectory.Text = Options.GenomeStarIndexDirectory;
            txtDbsnpVcfReference.Text = Options.ReferenceVcf;
            txtProteinFasta.Text = Options.ProteinFastaPath;
            txtStarFusionReference.Text = Options.Reference;
            txtUniProtProteinXml.Text = Options.UniProtXml;
        }

        private void PopulateChoices()
        {
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
                CbxWorkFlowType.Items.Add(aWorkFlow);
        }
    }
}