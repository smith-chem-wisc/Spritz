using CMD;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using ToolWrapperLayer;
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

            int ii = CmbxIndelFinding.SelectedIndex;
            if (ii == 0)
                Options.IndelFinder = "none";
            if (ii == 1)
                Options.IndelFinder = "gatk";
            if (ii == 2)
                Options.IndelFinder = "scalpel";

            //Options.SpritzDirectory = txtSpritzDirecory.Text;
            Options.AnalysisDirectory = txtAnalysisDirectory.Text;

            if (!Directory.Exists(Options.AnalysisDirectory))
            {
                MessageBox.Show("Analysis directory does not exist.", "Workflow", MessageBoxButton.OK);
                return;
            }

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
            Options.Reference = txtEnsemblReference.Text;
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
                else if (options.Command == LncRNADiscoveryFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 1;
                }
                else if (options.Command == TranscriptQuantificationFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 2;
                }
                else if (options.Command == STARAlignmentFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 3;
                }
                else if (options.Command == GeneFusionDiscoveryFlow.Command)
                {
                    CbxWorkFlowType.SelectedIndex = 4;
                }
                else
                {
                    // do nothing
                }
            }

            //txtSpritzDirecory.Text = options.SpritzDirectory;
            txtAnalysisDirectory.Text = AnalysisDirectory;
            txtThreads.Text = options.Threads.ToString();
            txtGenomeDir.Text = TrimQuotesOrNull(options.GenomeStarIndexDirectory);
            txtNewGeneModelGtfOrGff.Text = TrimQuotesOrNull(options.NewGeneModelGtfOrGff);
            txtEnsemblReference.Text = options.Reference ?? "GRCh38";
            txtUniProtProteinXml.Text = TrimQuotesOrNull(options.UniProtXml);
            ckbOverWriteStarAlignment.IsChecked = options.OverwriteStarAlignments;
            ckbStrandSpecific.IsChecked = options.StrandSpecific;
            ckbInferStrandedness.IsChecked = options.InferStrandSpecificity;
            CkbDoTranscriptIsoformAnalysis.IsChecked = options.DoTranscriptIsoformAnalysis;
            CkbDoGeneFusionAnalysis.IsChecked = options.DoFusionAnalysis;
            txtProteinFasta.Text = options.ProteinFastaPath;
            UpdateReference();
        }

        private string TrimQuotesOrNull(string a)
        {
            return a == null ? a : a.Trim('"');
        }

        private void PopulateChoices()
        {
            foreach (string aWorkFlow in Enum.GetNames(typeof(MyWorkflow)))
            {
                CbxWorkFlowType.Items.Add(aWorkFlow);
            }
            CmbxIndelFinding.Items.Add("None");
            CmbxIndelFinding.Items.Add("GATK");
            CmbxIndelFinding.Items.Add("Scalpel");
            CmbxIndelFinding.SelectedIndex = 2;
        }

        private void txtStarFusionReference_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            UpdateReference();
        }

        private void UpdateReference()
        {
            // Does dbSNP vcf already exist?
            var gatk = new GATKWrapper();
            if (gatk.KnownVariantSitesFileExists(EverythingRunnerEngine.SpritzDirectory, true, txtEnsemblReference.Text))
            {
                string ensemblVcfPath = Path.Combine(Path.GetDirectoryName(gatk.UcscKnownSitesPath), Path.GetFileNameWithoutExtension(gatk.UcscKnownSitesPath) + ".ensembl.vcf");
                txtDbsnpVcfReference.Text = ensemblVcfPath;
            }
            else
            {
                txtDbsnpVcfReference.Text = TrimQuotesOrNull(null);
            }

            // Does gene model already exist?
            if (AnalysisDirectory != null && AnalysisDirectory != "" && !Directory.Exists(AnalysisDirectory))
            {
                MessageBox.Show("Analysis directory does not exist.", "Workflow", MessageBoxButton.OK);
                return;
            }
            var ensembl = new EnsemblDownloadsWrapper();
            ensembl.DownloadReferences(EverythingRunnerEngine.SpritzDirectory, EverythingRunnerEngine.SpritzDirectory, txtEnsemblReference.Text, true);
            if (File.Exists(ensembl.Gff3GeneModelPath))
            {
                txtGeneModelGtfOrGff.Text = ensembl.Gff3GeneModelPath;
            }
            else if (File.Exists(ensembl.GtfGeneModelPath))
            {
                txtGeneModelGtfOrGff.Text = ensembl.GtfGeneModelPath;
            }
            else
            {
                txtGeneModelGtfOrGff.Text = TrimQuotesOrNull(null);
            }

            // Does genome reference already exist?
            if (File.Exists(ensembl.GenomeFastaPath))
            {
                txtGenomeFasta.Text = ensembl.GenomeFastaPath;
            }
            else
            {
                txtGenomeFasta.Text = TrimQuotesOrNull(null);
            }
        }
    }
}