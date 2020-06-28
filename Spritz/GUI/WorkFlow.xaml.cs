using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for workflows
    /// </summary>
    public partial class WorkFlowWindow : Window
    {
        private string AnalysisDirectory { get; set; }

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

        public Options Options { get; set; } = new Options();

        private MainWindow MainWindow { get; set; }

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
            Ensembl ensembl = (Ensembl)EnsemblReleaseVersions.SelectedItem;
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
            txtThreads.Text = options.Threads.ToString();
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

            EnsemblReleases = new ObservableCollection<Ensembl>();

            // read release.txt files into a list
            string releasefolder = Path.Combine(Directory.GetCurrentDirectory(), "EnsemblReleases");
            var releases = Directory.GetFiles(releasefolder, "*.txt").Select(Path.GetFileNameWithoutExtension).ToList();

            var genomeDB = File.ReadAllLines(Path.Combine(Directory.GetCurrentDirectory(), "genomes.csv"));
            foreach (string release in releases)
            {
                // read txt file into obsv collection
                var file = File.ReadAllLines(Path.Combine(releasefolder, $"{release}.txt"));
                var species = new List<string>(file);
                Dictionary<string, string> genomes = new Dictionary<string, string>();
                Dictionary<string, string> organisms = new Dictionary<string, string>();

                var unsupported = new List<string>();
                foreach (string genome in genomeDB.Where(g => g.Contains(release)))
                {
                    var splt = genome.Split(',');
                    genomes.Add(splt[1], splt[3]); // <Species, GenomeVer>
                    organisms.Add(splt[1], splt[2]); // <Species, OrganismName>

                    if (!string.Equals(splt[4], "86")) // only add species supported in snpeff (ver 86 ensembl)
                    {
                        unsupported.Add(splt[1]);
                    }
                }

                var supported = species.Where(s => !unsupported.Contains(s)).ToList();
                if (genomes.Count > 0)
                {
                    EnsemblReleases.Add(new Ensembl() { Release = release, Species = new ObservableCollection<string>(supported), Genomes = genomes, Organisms = organisms });
                }
            }
        }

        public ObservableCollection<Ensembl> EnsemblReleases { get; set; }

        public class Ensembl
        {
            public string Release { get; set; }
            public ObservableCollection<string> Species { get; set; }
            public Dictionary<string, string> Genomes { get; set; } // Mus_musculus GRCm38
            public Dictionary<string, string> Organisms { get; set; } // Mus_musculus GRCm38
        }

        private void txtThreads_TextChanged(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
        }

        public string Reference { get; set; } // define notify property changed

        private void Species_SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            saveButton.IsEnabled = true;

            // get selection from species
            var selectedEnsembl = (Ensembl)EnsemblReleaseVersions.SelectedItem;
            var selectedSpecies = (string)EnsemblSpecies.SelectedItem;
            Reference = selectedEnsembl.Genomes[selectedSpecies];
        }
    }
}