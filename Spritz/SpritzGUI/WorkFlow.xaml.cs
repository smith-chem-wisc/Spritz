using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using SpritzBackend;

namespace Spritz
{
    /// <summary>
    /// Interaction logic for workflows
    /// </summary>
    public partial class WorkFlowWindow : Window
    {
        private string AnalysisDirectory { get; set; }
        public string Reference { get; set; } // define notify property changed
        public ObservableCollection<EnsemblRelease> EnsemblReleases { get; set; }
        public SpritzOptions Options { get; set; } = new();
        private MainWindow MainWindow { get; set; }
        private int Threads { get; set; }

        public WorkFlowWindow(string analysisDirectory)
        {
            AnalysisDirectory = analysisDirectory;
            Options.AnalysisDirectory = SpritzOptions.DefaultAnalysisDirectory();
            Options.Threads = Environment.ProcessorCount;
            InitializeComponent();
            PopulateChoices();
            MainWindow = (MainWindow)Application.Current.MainWindow;
            UpdateFieldsFromTask(Options);
            DataContext = this;
        }

        public WorkFlowWindow(SpritzOptions options)
        {
            InitializeComponent();
            PopulateChoices();
            UpdateFieldsFromTask(options);
            MainWindow = (MainWindow)Application.Current.MainWindow;
            Options = options;
            DataContext = this;
        }

        protected void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        protected void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            //// Experiment type selection
            //int iii = CmbxExperimentType.SelectedIndex;
            //if (iii == 0)
            //{
            //    Options.ExperimentType = ExperimentType.RNASequencing.ToString();
            //}
            //else if (iii == 1)
            //{
            //    Options.ExperimentType = ExperimentType.WholeGenomeSequencing.ToString();
            //}
            //else if (iii == 2)
            //{
            //    Options.ExperimentType = ExperimentType.ExomeSequencing.ToString();
            //}
            //else
            //{
            //    MessageBox.Show("Please choose an experiment type selection.");
            //    return;
            //}

            Options.AnalysisDirectory = RunnerEngine.TrimQuotesOrNull(txtAnalysisDirectory.Text);
            if (!RunnerEngine.IsDirectoryWritable(Options.AnalysisDirectory))
            {
                MessageBox.Show($"Error: Cannot write to specified analysis directory: {Options.AnalysisDirectory}. Please choose another directory.",
                    "Write Permissions", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }
            Directory.CreateDirectory(Options.AnalysisDirectory);

            Options.Threads = Threads;
            EnsemblRelease ensembl = (EnsemblRelease)EnsemblReleaseVersions.SelectedItem;
            string species = EnsemblSpecies.SelectedItem.ToString();
            Options.Reference = EnsemblRelease.GetReferenceString(ensembl.Release, species, ensembl.Organisms[species], ensembl.Genomes[species]);
            Options.AnalyzeVariants = (bool)Cb_AnalyzeVariants.IsChecked;
            Options.AnalyzeIsoforms = (bool)Cb_AnalyzeIsoforms.IsChecked;
            Options.Quantify = (bool)Cb_Quantify.IsChecked;
            DialogResult = true;
        }

        private void UpdateFieldsFromTask(SpritzOptions options)
        {
            // Get information about the fastq and sra selections
            var rnaSeqFastqCollection = (ObservableCollection<RNASeqFastqDataGrid>)MainWindow.DataGridRnaSeqFastq.DataContext;
            Options.Fastq1 = string.Join(",", rnaSeqFastqCollection.Where(p => p.IsPairedEnd && p.MatePair == 1.ToString()).OrderBy(p => p.FileName).Select(p => p.FileName.Substring(0, p.FileName.Length - 2)).ToArray());
            Options.Fastq2 = string.Join(",", rnaSeqFastqCollection.Where(p => p.IsPairedEnd && p.MatePair == 2.ToString()).OrderBy(p => p.FileName).Select(p => p.FileName.Substring(0, p.FileName.Length - 2)).ToArray());
            Options.Fastq1SingleEnd = string.Join(",", rnaSeqFastqCollection.Where(p => !p.IsPairedEnd && p.MatePair == 1.ToString()).OrderBy(p => p.FileName).Select(p => p.FileName.Substring(0, p.FileName.Length - 2)).ToArray());

            var fq1s = Options.Fastq1.Split(',') ?? Array.Empty<string>();
            var fq2s = Options.Fastq2.Split(',') ?? Array.Empty<string>();
            var fq1s_se = Options.Fastq1SingleEnd.Split(',') ?? Array.Empty<string>();

            HashSet<string> unpairedFqPrefixes = new(
                fq1s.Where(fq1 => !fq2s.Any(fq2 => fq2.CompareTo(fq1) == 0)).Concat(
                    fq2s.Where(fq2 => !fq1s.Any(fq1s => fq1s.CompareTo(fq2) == 0))));
            if (unpairedFqPrefixes.Count > 0)
            {
                MessageBox.Show($"Add both paired files for {string.Join(",", unpairedFqPrefixes)}.",
                    "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Warning);
                throw new InvalidOperationException();
            }
            //HashSet<string> directories = new HashSet<string>(fq1s.Concat(fq2s).Concat(fq1s_se).Select(fq => Path.GetDirectoryName(fq)));
            //if (directories.Count > 1)
            //{
            //    MessageBox.Show($"All user-specified FASTQs must reside in one directory. Currently, they reside in {string.Join(",", directories)}.",
            //        "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Warning);
            //    throw new InvalidOperationException();
            //}

            //Options.ExperimentType = CmbxExperimentType.SelectedItem.ToString();
            var sraCollection = (ObservableCollection<SRADataGrid>)MainWindow.LbxSRAs.ItemsSource;
            Options.SraAccession = string.Join(",", sraCollection.Where(p => p.IsPairedEnd).Select(p => p.Name).ToArray());
            Options.SraAccessionSingleEnd = string.Join(",", sraCollection.Where(p => !p.IsPairedEnd).Select(p => p.Name).ToArray());
            if (Options.SraAccession.Length == 0 && options.Fastq1.Length == 0 && 
                Options.SraAccessionSingleEnd.Length == 0 && options.Fastq1SingleEnd.Length == 0)
            {
                Cb_AnalyzeIsoforms.IsChecked = false;
                Cb_AnalyzeIsoforms.IsEnabled = false;
                Cb_AnalyzeVariants.IsChecked = false;
                Cb_AnalyzeVariants.IsEnabled = false;
                Cb_Quantify.IsChecked = false;
                Cb_Quantify.IsEnabled = false;
            }

            txtAnalysisDirectory.Text = AnalysisDirectory;
            txtThreads.Text = MainWindow.DockerCPUs.ToString();
            Threads = MainWindow.DockerCPUs;
            Lb_ThreadInfo.Content = $"Integer between 1 and {MainWindow.DockerCPUs};\nmaximum is set in Docker Desktop";
            saveButton.IsEnabled = false;
        }

        private void PopulateChoices()
        {
            //CmbxExperimentType.Items.Add(ExperimentType.RNASequencing.ToString());
            //CmbxExperimentType.Items.Add(ExperimentType.WholeGenomeSequencing.ToString());
            //CmbxExperimentType.Items.Add(ExperimentType.ExomeSequencing.ToString());
            //CmbxExperimentType.SelectedIndex = 0; // hard coded selection (for now)

            EnsemblReleases = EnsemblRelease.GetReleases();
        }

        private void Species_SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            saveButton.IsEnabled = true;

            // get selection from species
            var selectedEnsembl = (EnsemblRelease)EnsemblReleaseVersions.SelectedItem;
            var selectedSpecies = (string)EnsemblSpecies.SelectedItem;
            Reference = selectedEnsembl.Genomes[selectedSpecies];
        }

        private void TxtThreads_LostFocus(object sender, RoutedEventArgs e)
        {
            if (int.TryParse(txtThreads.Text, out int threads) && threads <= MainWindow.DockerCPUs && threads > 0)
            {
                Threads = threads;
            }
            else
            {
                txtThreads.Text = MainWindow.DockerCPUs.ToString();
            }
        }
    }
}