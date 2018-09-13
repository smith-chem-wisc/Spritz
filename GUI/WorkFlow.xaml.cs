using CMD;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using WorkflowLayer;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for workflows
    /// </summary>
    public partial class WorkFlowWindow : Window
    {
        private string AnalysisDirectory;

        public WorkFlowWindow(string analysisDirectory)
        {
            AnalysisDirectory = analysisDirectory;
            InitializeComponent();
            PopulateChoices();
            Options = new Options();
            MainWindow = (MainWindow)Application.Current.MainWindow;
            UpdateFieldsFromTask(Options);
        }

        public WorkFlowWindow(Options options)
        {
            InitializeComponent();
            PopulateChoices();
            UpdateFieldsFromTask(options);
            MainWindow = (MainWindow)Application.Current.MainWindow;
            Options = options;
        }

        public Options Options { get; set; }

        private MainWindow MainWindow { get; set; }

        protected void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        protected void saveButton_Click(object sender, RoutedEventArgs e)
        {
            int i = CbxWorkFlowType.SelectedIndex;
            if (i == 0)
                Options.Command = SampleSpecificProteinDBFlow.Command;
            else if (i == 1)
                Options.Command = LncRNADiscoveryFlow.Command;
            else if (i == 2)
                Options.Command = TranscriptQuantificationFlow.Command;
            else if (i == 3)
                Options.Command = STARAlignmentFlow.Command;
            else if (i == 4)
                Options.Command = GeneFusionDiscoveryFlow.Command;
            else
            {
                MessageBox.Show("Please choose a workflow.");
                return;
            }

            //Options.SpritzDirectory = txtSpritzDirecory.Text;
            Options.AnalysisDirectory = txtAnalysisDirectory.Text;

            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)MainWindow.dataGridRnaSeqFastq.DataContext;
            if (rnaSeqFastqCollection.Count != 0)
            {
                Options.Fastq1 = string.Join(",", rnaSeqFastqCollection.Where(p => p.MatePair == 1.ToString()).OrderBy(p => p.Experiment).Select(p => p.FilePath).ToArray());
                Options.Fastq2 = string.Join(",", rnaSeqFastqCollection.Where(p => p.MatePair == 2.ToString()).OrderBy(p => p.Experiment).Select(p => p.FilePath).ToArray());
            }
            var sraCollection = (ObservableCollection<SRADataGrid>)MainWindow.LbxSRAs.ItemsSource;
            Options.SraAccession = string.Join(",", sraCollection.Select(p => p.Name).ToArray());

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

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (Directory.Exists(draggedFilePath))
                    {
                        foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                        {
                            AddAFile(file);
                        }
                    }
                    else
                    {
                        AddAFile(draggedFilePath);
                    }
                }
            }
        }

        private void AddAFile(string filepath)
        {
            var theExtension = filepath.EndsWith("gz") ?
                Path.GetExtension(Path.GetExtension(filepath)).ToLowerInvariant() :
                Path.GetExtension(filepath).ToLowerInvariant();

            if (theExtension == ".fa")
                txtGenomeFasta.Text = filepath;
            else if (theExtension == ".xml")
                txtUniProtProteinXml.Text = filepath;
            else if (theExtension == ".gtf" || theExtension.Contains("gff"))
                txtGeneModelGtfOrGff.Text = filepath;
            else if (theExtension == ".vcf")
                txtDbsnpVcfReference.Text = filepath;
            else
                return;
        }

        private void UpdateFieldsFromTask(Options options)
        {
            foreach (var aWorkFlow in Enum.GetValues(typeof(MyWorkflow)))
            {
                if (options.Command == SampleSpecificProteinDBFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 0;
                }
                if (options.Command == LncRNADiscoveryFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 1;
                }
                if (options.Command == TranscriptQuantificationFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 2;
                }
                if (options.Command == STARAlignmentFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 3;
                }
                if (options.Command == GeneFusionDiscoveryFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 4;
                }
            }

            //txtSpritzDirecory.Text = options.SpritzDirectory;
            txtAnalysisDirectory.Text = AnalysisDirectory;

            txtThreads.Text = options.Threads.ToString();
            txtGenomeDir.Text = TrimQuotesOrNull(options.GenomeStarIndexDirectory);
            txtGenomeFasta.Text = TrimQuotesOrNull(options.GenomeFasta);
            txtGeneModelGtfOrGff.Text = TrimQuotesOrNull(options.GeneModelGtfOrGff);
            txtNewGeneModelGtfOrGff.Text = TrimQuotesOrNull(options.NewGeneModelGtfOrGff);
            txtDbsnpVcfReference.Text = TrimQuotesOrNull(options.ReferenceVcf);
            txtStarFusionReference.Text = options.Reference;
            txtUniProtProteinXml.Text = TrimQuotesOrNull(options.UniProtXml);
            ckbOverWriteStarAlignment.IsChecked = options.OverwriteStarAlignments;
            ckbStrandSpecific.IsChecked = options.StrandSpecific;
            ckbInferStrandedness.IsChecked = options.InferStrandSpecificity;
            CkbDoTranscriptIsoformAnalysis.IsChecked = options.DoTranscriptIsoformAnalysis;
            CkbDoGeneFusionAnalysis.IsChecked = options.DoFusionAnalysis;
            txtProteinFasta.Text = options.ProteinFastaPath;
        }

        private string TrimQuotesOrNull(string a)
        {
            return a == null ? a : a.Trim('"');
        }

        private void PopulateChoices()
        {
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
                CbxWorkFlowType.Items.Add(aWorkFlow);
        }
    }
}