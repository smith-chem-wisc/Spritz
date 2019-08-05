using CMD;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using ToolWrapperLayer;
using WorkflowLayer;
using PerformanceCounter = System.Diagnostics.PerformanceCounter;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for workflows
    /// </summary>
    public partial class WorkFlowWindow : Window
    {
        private string AnalysisDirectory { get; set; }
        private float AvailableMemoryMb { get; set; }
        private bool SettingUp { get; set; } = true;
        private static readonly int GatkMemoryReccommendation = 16000;
        private static readonly int ScalpelMemoryReccommendation = 32000;

        public WorkFlowWindow(string analysisDirectory)
        {
            AnalysisDirectory = analysisDirectory;
            InitializeComponent();
            PopulateChoices();
            MainWindow = (MainWindow)Application.Current.MainWindow;
            UpdateFieldsFromTask(Options);
            SettingUp = false;
            ChooseWorkerPreset();
        }

        public WorkFlowWindow(Options options)
        {
            InitializeComponent();
            PopulateChoices();
            UpdateFieldsFromTask(options);
            MainWindow = (MainWindow)Application.Current.MainWindow;
            Options = options;
            SettingUp = false;
            ChooseWorkerPreset();
        }

        public Options Options { get; set; } = new Options();

        private MainWindow MainWindow { get; set; }

        protected void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        protected void SaveButton_Click(object sender, RoutedEventArgs e)
        {           
            // Command selection
            int i = CbxWorkFlowType.SelectedIndex;
            if (i == 0)
            {
                Options.Command = SampleSpecificProteinDBFlow.Command;
            }
            else if (i == 1)
            {
                Options.Command = LncRNADiscoveryFlow.Command;
            }
            else if (i == 2)
            {
                Options.Command = TranscriptQuantificationFlow.Command;
            }
            else if (i == 3)
            {
                Options.Command = AlignmentFlow.Command;
            }
            else if (i == 4)
            {
                Options.Command = GeneFusionDiscoveryFlow.Command;
            }
            else
            {
                MessageBox.Show("Please choose a workflow.");
                return;
            }

            // Indel finder selection
            int ii = CmbxIndelFinding.SelectedIndex;
            if (ii == 0)
            {
                Options.IndelFinder = "none";
            }
            else if (ii == 1)
            {
                Options.IndelFinder = "gatk";
            }
            else if (ii == 2)
            {
                Options.IndelFinder = "scalpel";
            }
            else
            {
                MessageBox.Show("Please choose an indel finding selection.");
                return;
            }

            // Experiment type selection
            int iii = CmbxExperimentType.SelectedIndex;
            if (iii == 0)
            {
                Options.ExperimentType = ExperimentType.RNASequencing.ToString();
            }
            else if (iii == 1)
            {
                Options.ExperimentType = ExperimentType.WholeGenomeSequencing.ToString();
            }
            else if (iii == 2)
            {
                Options.ExperimentType = ExperimentType.ExomeSequencing.ToString();
            }
            else
            {
                MessageBox.Show("Please choose an experiment type selection.");
                return;
            }

            // Options.SpritzDirectory = txtSpritzDirecory.Text;
            var defaultAnalysisDirectory = Path.Combine(Directory.GetCurrentDirectory(), "output");
            if (!Directory.Exists(defaultAnalysisDirectory))
            {
                Directory.CreateDirectory(defaultAnalysisDirectory);
            }

            Options.AnalysisDirectory = TrimQuotesOrNull(txtAnalysisDirectory.Text);

            if (!Directory.Exists(Options.AnalysisDirectory))
            {

                MessageBox.Show("Analysis directory does not exist.", "Workflow", MessageBoxButton.OK);
                return;
            }

            Options.Threads = int.Parse(txtThreads.Text);

            // features yet to be supported
            //Options.GenomeStarIndexDirectory = txtGenomeDir.Text;
            //Options.GenomeFasta = txtGenomeFasta.Text;
            //Options.GeneModelGtfOrGff = txtGeneModelGtfOrGff.Text;
            //Options.NewGeneModelGtfOrGff = txtNewGeneModelGtfOrGff.Text;
            //Options.ReferenceVcf = txtDbsnpVcfReference.Text;
            //Options.Reference = txtEnsemblReference.Text;
            //Options.UniProtXml = txtUniProtProteinXml.Text;
            //Options.OverwriteStarAlignments = ckbOverWriteStarAlignment.IsChecked.Value;
            //Options.StrandSpecific = ckbStrandSpecific.IsChecked.Value;
            //Options.InferStrandSpecificity = ckbInferStrandedness.IsChecked.Value;
            //Options.SkipVariantAnalysis = CkbSkipVariantAnalysis.IsChecked.Value;
            //Options.DoTranscriptIsoformAnalysis = CkbDoTranscriptIsoformAnalysis.IsChecked.Value;
            //Options.DoFusionAnalysis = CkbDoGeneFusionAnalysis.IsChecked.Value;
            //Options.VariantCallingWorkers = int.Parse(TxtVariantCallingWorkerNum.Text);
            //Options.ProteinFastaPath = txtProteinFasta.Text;
            DialogResult = true;
        }

        //private void Window_Drop(object sender, DragEventArgs e)
        //{
        //    string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
        //    if (files != null)
        //    {
        //        foreach (var draggedFilePath in files)
        //        {
        //            if (Directory.Exists(draggedFilePath))
        //            {
        //                foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
        //                {
        //                    AddAFile(file);
        //                }
        //            }
        //            else
        //            {
        //                AddAFile(draggedFilePath);
        //            }
        //        }
        //    }
        // }

        //private void AddAFile(string filepath)
        //{
        //    var theExtension = filepath.EndsWith("gz") ?
        //        Path.GetExtension(Path.GetFileNameWithoutExtension(filepath)).ToLowerInvariant() :
        //        Path.GetExtension(filepath).ToLowerInvariant();

        //    if (theExtension == ".fa")
        //        txtGenomeFasta.Text = filepath;
        //    else if (theExtension == ".xml")
        //        txtUniProtProteinXml.Text = filepath;
        //    else if (theExtension == ".gtf" || theExtension.Contains("gff"))
        //        txtGeneModelGtfOrGff.Text = filepath;
        //    else if (theExtension == ".vcf")
        //        txtDbsnpVcfReference.Text = filepath;
        //    else
        //        return;
        //}

        private void UpdateFieldsFromTask(Options options)
        {
            // Get information about the fastq and sra selections
            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)MainWindow.DataGridRnaSeqFastq.DataContext;
            Options.Fastq1 = string.Join(",", rnaSeqFastqCollection.Where(p => p.MatePair == 1.ToString()).OrderBy(p => p.FileName).Select(p => p.FileName.Substring(0, p.FileName.Length - 2)).ToArray());
            Options.Fastq2 = string.Join(",", rnaSeqFastqCollection.Where(p => p.MatePair == 2.ToString()).OrderBy(p => p.FileName).Select(p => p.FileName.Substring(0, p.FileName.Length - 2)).ToArray());

            //use RNAFastqCollection instead
            var fq1s = Options.Fastq1.Split(',') ?? new string[0];
            var fq2s = Options.Fastq2.Split(',') ?? new string[0];

            foreach (string fq1 in fq1s)
            {
                if (!fq2s.Any(fq2 => fq2.CompareTo(fq1) == 0))
                {
                    MessageBox.Show("Only paired end sequencing is supported. Add both paired files for " + fq1 + ".", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                    throw new InvalidOperationException();
                }
            }

            Options.ExperimentType = CmbxExperimentType.SelectedItem.ToString();
            var sraCollection = (ObservableCollection<SRADataGrid>)MainWindow.LbxSRAs.ItemsSource;
            Options.SraAccession = string.Join(",", sraCollection.Select(p => p.Name).ToArray());

            // add the commands
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
                else if (options.Command == AlignmentFlow.Command)
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

            txtAnalysisDirectory.Text = AnalysisDirectory;
            txtThreads.Text = options.Threads.ToString();
            txtEnsemblReference.Text = options.Reference ?? "GRCh38";
            TxtVariantCallingWorkerNum.Text = options.VariantCallingWorkers.ToString();

            //txtSpritzDirecory.Text = options.SpritzDirectory;
            //txtGenomeDir.Text = TrimQuotesOrNull(options.GenomeStarIndexDirectory);
            //txtNewGeneModelGtfOrGff.Text = TrimQuotesOrNull(options.NewGeneModelGtfOrGff);
            //txtUniProtProteinXml.Text = TrimQuotesOrNull(options.UniProtXml);
            //ckbOverWriteStarAlignment.IsChecked = options.OverwriteStarAlignments;
            //ckbStrandSpecific.IsChecked = options.StrandSpecific;
            //ckbInferStrandedness.IsChecked = options.InferStrandSpecificity;
            //CkbSkipVariantAnalysis.IsChecked = options.Fastq1 == null && options.SraAccession == null || options.SkipVariantAnalysis;
            //CkbDoTranscriptIsoformAnalysis.IsChecked = options.DoTranscriptIsoformAnalysis;
            //CkbDoGeneFusionAnalysis.IsChecked = options.DoFusionAnalysis;
            //txtProteinFasta.Text = options.ProteinFastaPath;
            //UpdateReference();
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
            if (Options.Fastq1 == null)
            {
                Options.Command = SampleSpecificProteinDBFlow.Command; // this is the only one currently allowed without FASTQ files
            }

            CmbxIndelFinding.Items.Add("None");
            CmbxIndelFinding.Items.Add("GATK");
            CmbxIndelFinding.Items.Add("Scalpel");
            CmbxIndelFinding.SelectedIndex = 1; // hard coded selection (for now)

            CmbxExperimentType.Items.Add(ExperimentType.RNASequencing.ToString());
            CmbxExperimentType.Items.Add(ExperimentType.WholeGenomeSequencing.ToString());
            CmbxExperimentType.Items.Add(ExperimentType.ExomeSequencing.ToString());
            CmbxExperimentType.SelectedIndex = 0; // hard coded selection (for now)
        }

        //private void txtStarFusionReference_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        //{
        //    UpdateReference();
        //}

        //private void UpdateReference()
        //{
        //    // Does dbSNP vcf already exist?
        //    var gatk = new GATKWrapper(1);
        //    string ensemblVcfPath = gatk.DownloadEnsemblKnownVariantSites(EverythingRunnerEngine.SpritzDirectory, true, txtEnsemblReference.Text, true);
        //    if (File.Exists(ensemblVcfPath))
        //    {
        //        txtDbsnpVcfReference.Text = ensemblVcfPath;
        //    }
        //    else
        //    {
        //        txtDbsnpVcfReference.Text = TrimQuotesOrNull(null);
        //    }

        //    // Does gene model already exist?
        //    if (AnalysisDirectory != null && AnalysisDirectory != "" && !Directory.Exists(AnalysisDirectory))
        //    {
        //        MessageBox.Show("Analysis directory does not exist.", "Workflow", MessageBoxButton.OK);
        //        return;
        //    }
        //    var ensembl = new EnsemblDownloadsWrapper();
        //    ensembl.DownloadReferences(EverythingRunnerEngine.SpritzDirectory, EverythingRunnerEngine.SpritzDirectory, txtEnsemblReference.Text, true);
        //    if (File.Exists(ensembl.Gff3GeneModelPath))
        //    {
        //        txtGeneModelGtfOrGff.Text = ensembl.Gff3GeneModelPath;
        //    }
        //    else if (File.Exists(ensembl.GtfGeneModelPath))
        //    {
        //        txtGeneModelGtfOrGff.Text = ensembl.GtfGeneModelPath;
        //    }
        //    else
        //    {
        //        txtGeneModelGtfOrGff.Text = TrimQuotesOrNull(null);
        //    }

        //    // Does genome reference already exist?
        //    if (File.Exists(ensembl.GenomeFastaPath))
        //    {
        //        txtGenomeFasta.Text = ensembl.GenomeFastaPath;
        //    }
        //    else
        //    {
        //        txtGenomeFasta.Text = TrimQuotesOrNull(null);
        //    }
        //}

        private void ChooseWorkerPreset()
        {
            if (SettingUp)
            {
                return;
            }

            if (AvailableMemoryMb <= 0)
            {
                MessageBox.Show("Choosing some presets. This will take a couple seconds. Hit OK.", "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
                var performance = new PerformanceCounter("Memory", "Available MBytes");
                AvailableMemoryMb = performance.NextValue();
            }

            string indelFinder = CmbxIndelFinding.Items[CmbxIndelFinding.SelectedIndex].ToString(); // get new selected index (.Text gives old result)
            int recommendation = "scalpel".Equals(indelFinder, StringComparison.InvariantCultureIgnoreCase) ? 32000 : 16000;
            bool validThreads = int.TryParse(txtThreads.Text, out int threads);

            int workers = 1;
            if (validThreads && threads == 1)
            {
                // do nothing if there's just one thread
            }
            else
            {
                float workerMemory = AvailableMemoryMb;
                int workerThreads = threads;
                while (workerThreads > 2 && workerMemory > recommendation)
                {
                    workers++;
                    workerMemory = AvailableMemoryMb / (float)workers;
                    workerThreads = (int)Math.Ceiling((float)threads / (float)workers);
                }
                workers--;
            }
            TxtVariantCallingWorkerNum.Text = workers.ToString();
        }

        private void CheckMemory()
        {
            if (SettingUp)
            {
                return;
            }

            bool validThreads = int.TryParse(txtThreads.Text, out int threads);
            bool validVariantCallers = int.TryParse(TxtVariantCallingWorkerNum.Text, out int variantCallingWorkers);
            if (validThreads
                && validVariantCallers
                && variantCallingWorkers > 0
                && variantCallingWorkers <= Math.Ceiling((double)threads / (double)2))
            {
                Options.VariantCallingWorkers = variantCallingWorkers;
            }

            string indelFinder = CmbxIndelFinding.Items[CmbxIndelFinding.SelectedIndex].ToString(); // get new selected index (.Text gives old result)
            int recommendation = "scalpel".Equals(indelFinder, StringComparison.InvariantCultureIgnoreCase) ? ScalpelMemoryReccommendation : GatkMemoryReccommendation;
            if (Options.Command == SampleSpecificProteinDBFlow.Command
                && AvailableMemoryMb / Options.VariantCallingWorkers < recommendation)
            {
                MessageBox.Show($"Using {Options.VariantCallingWorkers.ToString()} workers with {AvailableMemoryMb.ToString()} MB of memory " +
                    $"leaves less than the {recommendation} MB memory per worker recommended for variant calling.",
                    "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
            }

            if (variantCallingWorkers > threads)
            {
                MessageBox.Show($"Using {Options.VariantCallingWorkers.ToString()} workers with {threads.ToString()} threads " +
                    $"leaves less than the {2} threads per worker recommended for variant calling.",
                    "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
            }

            if (validThreads
                && validVariantCallers
                && variantCallingWorkers <= 0)
            {
                TxtVariantCallingWorkerNum.Text = (++variantCallingWorkers).ToString(); 
            }
        }

        private void CmdUpVariantCallingWorkers_Click(object sender, RoutedEventArgs e)
        {
            TxtVariantCallingWorkerNum.Text = (++Options.VariantCallingWorkers).ToString();
        }

        private void CmdDownVariantCallingWorkers_Click(object sender, RoutedEventArgs e)
        {
            TxtVariantCallingWorkerNum.Text = (--Options.VariantCallingWorkers).ToString();
        }

        private void TxtVariantCallingWorkerNum_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            CheckMemory();
        }

        private void CmbxIndelFinding_SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            ChooseWorkerPreset();
        }

        private void txtThreads_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            ChooseWorkerPreset();
        }
    }
}