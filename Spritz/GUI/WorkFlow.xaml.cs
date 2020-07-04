using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Security;
using System.Security.Permissions;
using System.Windows;

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
        public Options Options { get; set; } = new Options(Environment.ProcessorCount);
        private MainWindow MainWindow { get; set; }
        private int Threads { get; set; }

        public WorkFlowWindow(string analysisDirectory)
        {
            AnalysisDirectory = analysisDirectory;
            InitializeComponent();
            PopulateChoices();
            MainWindow = (MainWindow)Application.Current.MainWindow;
            UpdateFieldsFromTask(Options);
            DataContext = this;
        }

        public WorkFlowWindow(Options options)
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

            Options.AnalysisDirectory = TrimQuotesOrNull(txtAnalysisDirectory.Text);
            try
            {
                string testDirectory = Path.Combine(Options.AnalysisDirectory, $"TestSpritzPermissions{Options.AnalysisDirectory.GetHashCode()}");
                Directory.CreateDirectory(testDirectory);
                Directory.Delete(testDirectory);
            }
            catch (Exception)
            {
                MessageBox.Show($"Error: Cannot write to specified analysis directory: {Options.AnalysisDirectory}. Please choose another directory.", "Write Permissions", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }

            if (!Directory.Exists(Options.AnalysisDirectory))
            {
                Directory.CreateDirectory(Options.AnalysisDirectory);
            }

            Options.Threads = Threads;
            EnsemblRelease ensembl = (EnsemblRelease)EnsemblReleaseVersions.SelectedItem;
            Options.Release = ensembl.Release;
            Options.Species = EnsemblSpecies.SelectedItem.ToString();
            Options.Reference = ensembl.Genomes[Options.Species];
            Options.Organism = ensembl.Organisms[Options.Species];
            Options.SnpEff = "86";
            DialogResult = true;
        }

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

            //Options.ExperimentType = CmbxExperimentType.SelectedItem.ToString();
            var sraCollection = (ObservableCollection<SRADataGrid>)MainWindow.LbxSRAs.ItemsSource;
            Options.SraAccession = string.Join(",", sraCollection.Select(p => p.Name).ToArray());
            txtAnalysisDirectory.Text = AnalysisDirectory;
            txtThreads.Text = MainWindow.DockerCPUs.ToString();
            Threads = MainWindow.DockerCPUs;
            Lb_ThreadInfo.Content = $"Integer between 1 and {MainWindow.DockerCPUs};\nmaximum is set in Docker Desktop";
            saveButton.IsEnabled = false;
        }

        private string TrimQuotesOrNull(string a)
        {
            return a == null ? a : a.Trim('"');
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

        private void txtThreads_LostFocus(object sender, RoutedEventArgs e)
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