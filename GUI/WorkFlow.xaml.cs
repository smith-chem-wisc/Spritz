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
            UpdateFieldsFromTask(Options);
        }

        public WorkFlowWindow(Options options)
        {
            InitializeComponent();
            PopulateChoices();
            UpdateFieldsFromTask(options);
            mainWindow = (MainWindow)Application.Current.MainWindow;
            Options = options;
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
                    Options.Command = STARAlignmentFlow.Command;
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
                Options.Fastq1 = String.Join(",", rnaSeqFastqCollection.Where(p => p.MateRun == 1.ToString()).OrderBy(p=>p.Experiment).Select(p => p.FilePath).ToArray());
                Options.Fastq2 = String.Join(",", rnaSeqFastqCollection.Where(p => p.MateRun == 2.ToString()).OrderBy(p => p.Experiment).Select(p => p.FilePath).ToArray());
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
            Options.QuickSnpEffWithoutStats = CkbQuickSnpEffWithoutStats.IsChecked.Value;
            Options.ProteinFastaPath = txtProteinFasta.Text;          
            DialogResult = true;
        }

        private void UpdateFieldsFromTask(Options options)
        {
            foreach (var aWorkFlow in Enum.GetValues(typeof(MyWorkflow)))
            {
                if (options.Command == "proteins")
                {
                    CbxWorkFlowType.SelectedIndex = 0;
                }
                if (options.Command == "lncRNADiscovery")
                {
                    CbxWorkFlowType.SelectedIndex = 1;
                }
                if (options.Command == "quantify")
                {
                    CbxWorkFlowType.SelectedIndex = 2;
                }
                if (options.Command == "a")
                {
                    CbxWorkFlowType.SelectedIndex = 3;
                }
                if (options.Command == "fusion")
                {
                    CbxWorkFlowType.SelectedIndex = 4;
                }
            }

            txtSpritzDirecory.Text = options.SpritzDirectory;
            txtAnalysisDirectory.Text = options.AnalysisDirectory;

            txtThreads.Text = options.Threads.ToString();
            txtGenomeDir.Text = options.GenomeStarIndexDirectory;
            txtGenomeFasta.Text = options.GenomeFasta;
            txtGeneModelGtfOrGff.Text = options.GeneModelGtfOrGff;
            txtNewGeneModelGtfOrGff.Text = options.NewGeneModelGtfOrGff;
            txtDbsnpVcfReference.Text = options.ReferenceVcf;
            txtStarFusionReference.Text = options.Reference;
            txtUniProtProteinXml.Text = options.UniProtXml;
            ckbOverWriteStarAlignment.IsChecked = options.OverwriteStarAlignments;
            ckbStrandSpecific.IsChecked = options.StrandSpecific;
            ckbInferStrandedness.IsChecked = options.InferStrandSpecificity;
            CkbDoTranscriptIsoformAnalysis.IsChecked = options.DoTranscriptIsoformAnalysis;
            CkbDoGeneFusionAnalysis.IsChecked = options.DoFusionAnalysis;
            CkbQuickSnpEffWithoutStats.IsChecked = options.QuickSnpEffWithoutStats;
            txtProteinFasta.Text = options.ProteinFastaPath;
        }

        private void PopulateChoices()
        {
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
                CbxWorkFlowType.Items.Add(aWorkFlow);
        }
    }
}