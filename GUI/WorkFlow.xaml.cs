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

            Options.SpritzDirectory = txtSpritzDirecory.Text;
            Options.AnalysisDirectory = txtAnalysisDirectory.Text;    
            
            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)mainWindow.dataGridRnaSeqFastq.DataContext;
            if (rnaSeqFastqCollection.Count != 0)
            {
                Options.Fastq1 = rnaSeqFastqCollection.Select(p => p.FilePath).ToList()[0];
                Options.Fastq2 = rnaSeqFastqCollection.Select(p => p.FilePath).ToList()[1];
            }            
            var sraCollection = (ObservableCollection<SRADataGrid>)mainWindow.LbxSRAs.ItemsSource;
            Options.SraAccession = String.Join(",", sraCollection.Select(p => p.Name).ToArray());

            Options.Threads = int.Parse(txtThreads.Text);
            Options.GenomeStarIndexDirectory = txtGenomeDir.Text;
            Options.GenomeFasta = txtGenomeFasta.Text;
            Options.GeneModelGtfOrGff = txtGeneModelGtfOrGff.Text;
            Options.NewGeneModelGtfOrGff = txtNewGeneModelGtfOrGff.Text;
            Options.ReferenceVcf = txtDbsnpVcfReference.Text;
            Options.Reference = txtStarFusionReference.Text;
            Options.UniProtXml = txtUniProtProteinXml.Text;
            Options.OverwriteStarAlignments = ckbOverWriteStarAlignment.IsChecked.Value;
            Options.StrandSpecific = ckbStrandSpecific.IsChecked.Value;
            Options.InferStrandSpecificity = ckbInferStrandedness.IsChecked.Value;
            Options.DoTranscriptIsoformAnalysis = CkbDoTranscriptIsoformAnalysis.IsChecked.Value;
            Options.DoFusionAnalysis = CkbDoGeneFusionAnalysis.IsChecked.Value;

            Options.ProteinFastaPath = txtProteinFasta.Text;          
            DialogResult = true;
        }

        private void UpdateFieldsFromTask()
        {
            int i = 0;
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
            {
                if (Options.Command == aWorkFlow)
                {
                    CbxWorkFlowType.SelectedIndex = i;
                }
                i++;
            }

            txtSpritzDirecory.Text = Options.SpritzDirectory;
            txtAnalysisDirectory.Text = Options.AnalysisDirectory;

            txtThreads.Text = Options.Threads.ToString();
            txtGenomeDir.Text = Options.GenomeStarIndexDirectory;
            txtGenomeFasta.Text = Options.GenomeFasta;
            txtGeneModelGtfOrGff.Text = Options.GeneModelGtfOrGff;
            txtNewGeneModelGtfOrGff.Text = Options.NewGeneModelGtfOrGff;
            txtDbsnpVcfReference.Text = Options.ReferenceVcf;
            txtStarFusionReference.Text = Options.Reference;
            txtUniProtProteinXml.Text = Options.UniProtXml;
            ckbOverWriteStarAlignment.IsChecked = Options.OverwriteStarAlignments;
            ckbStrandSpecific.IsChecked = Options.StrandSpecific;
            ckbInferStrandedness.IsChecked = Options.InferStrandSpecificity;
            CkbDoTranscriptIsoformAnalysis.IsChecked = Options.DoTranscriptIsoformAnalysis;
            CkbDoGeneFusionAnalysis.IsChecked = Options.DoFusionAnalysis;

            txtProteinFasta.Text = Options.ProteinFastaPath;
        }

        private void PopulateChoices()
        {
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
                CbxWorkFlowType.Items.Add(aWorkFlow);
        }
    }
}