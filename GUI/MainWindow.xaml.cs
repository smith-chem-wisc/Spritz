using CMD;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using WorkflowLayer;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        #region Private Fields

        private readonly ObservableCollection<GenomeFastaDataGrid> genomeFastaCollection = new ObservableCollection<GenomeFastaDataGrid>();
        private readonly ObservableCollection<GeneSetDataGrid> geneSetCollection = new ObservableCollection<GeneSetDataGrid>();
        private readonly ObservableCollection<RNASeqFastqDataGrid> rnaSeqFastqCollection = new ObservableCollection<RNASeqFastqDataGrid>();
        private ObservableCollection<InRunTask> dynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
        private readonly ObservableCollection<PreRunTask> staticTasksObservableCollection = new ObservableCollection<PreRunTask>();

        #endregion Private Fields

        #region Public Constructors

        public MainWindow()
        {
            InitializeComponent();

            dataGridFASTA.DataContext = genomeFastaCollection;
            dataGridGeneSet.DataContext = geneSetCollection;
            dataGridRnaSeqFastq.DataContext = rnaSeqFastqCollection;
            workflowTreeView.DataContext = staticTasksObservableCollection;
        }

        #endregion Public Constructors

        #region Private Methods - Controlers

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            if (files != null)
                foreach (var draggedFilePath in files)
                {
                    if (Directory.Exists(draggedFilePath))
                        foreach (string file in Directory.EnumerateFiles(draggedFilePath, "*.*", SearchOption.AllDirectories))
                        {
                            AddAFile(file);
                        }
                    else
                    {
                        AddAFile(draggedFilePath);
                    }
                    dataGridFASTA.Items.Refresh();
                    dataGridGeneSet.Items.Refresh();
                    dataGridRnaSeqFastq.Items.Refresh();
                }
            UpdateOutputFolderTextbox();
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
        }

        private void MenuItem_Wiki_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/Spritz/wiki");
        }

        private void MenuItem_Contact_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"https://github.com/smith-chem-wisc/Spritz");
        }

        private void BtnSTARAlignment_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new TransferModificationsFlowWindows();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.Options);
                UpdateTaskGuiStuff();
            }
        }

        private void BtnAddFastq2Proteins_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new Fastq2ProteinsFlowWindows();
            if (dialog.ShowDialog() == true)
            {
                //AddTaskToCollection(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
        }

        private void BtnAddLncRNADiscover_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new TransferModificationsFlowWindows();
            if (dialog.ShowDialog() == true)
            {
                AddTaskToCollection(dialog.Options);
                UpdateTaskGuiStuff();
            }
        }

        private void RunWorkflowButton_Click(object sender, RoutedEventArgs e)
        {
            dynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
            for (int i = 0; i < staticTasksObservableCollection.Count; i++)
            {
                dynamicTasksObservableCollection.Add(new InRunTask("Workflow" + (i + 1) + "-" + staticTasksObservableCollection[i].options.Command.ToString(), staticTasksObservableCollection[i].options));
            }
            workflowTreeView.DataContext = dynamicTasksObservableCollection;
            EverythingRunnerEngine a = new EverythingRunnerEngine(dynamicTasksObservableCollection.Select(b => new Tuple<string, Options>(b.DisplayName, b.options)).ToList(), OutputFolderTextBox.Text);
            a.Run();
            //var t = new Task(a.Run);
            //t.Start();
        }

        private void DataGridCell_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void BtnAddGenomeFASTA_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Genome Fasta Files (*.fa; *.fasta)|*.fa; *.fasta",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var filepath in openPicker.FileNames)
                {
                    AddAFile(filepath);
                }
            dataGridFASTA.Items.Refresh();
        }

        private void BtnClearGenomeFASTA_Click(object sender, RoutedEventArgs e)
        {
            genomeFastaCollection.Clear();
        }

        private void BtnAddRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "FASTQ Files|*.fastq",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var filepath in openPicker.FileNames)
                {
                    AddAFile(filepath);
                }
            dataGridRnaSeqFastq.Items.Refresh();
        }

        private void BtnClearRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            rnaSeqFastqCollection.Clear();
        }

        private void BtnAddGeneSet_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Gene Model Files|*.gtf;*.gff3",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var filepath in openPicker.FileNames)
                {
                    AddAFile(filepath);
                }
            dataGridGeneSet.Items.Refresh();
        }

        private void BtnClearGeneSet_Click(object sender, RoutedEventArgs e)
        {
            geneSetCollection.Clear();
        }

        private void LoadTaskButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "TOML files(*.toml)|*.toml",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = true
            };
            if (openPicker.ShowDialog() == true)
                foreach (var tomlFromSelected in openPicker.FileNames)
                {
                    AddAFile(tomlFromSelected);
                }
            UpdateTaskGuiStuff();
        }

        private void ClearTasksButton_Click(object sender, RoutedEventArgs e)
        {
            staticTasksObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void RemoveLastTaskButton_Click(object sender, RoutedEventArgs e)
        {
            staticTasksObservableCollection.RemoveAt(staticTasksObservableCollection.Count - 1);
            UpdateTaskGuiStuff();
        }

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
        }

        private void AddNewRnaSeqFastq(object sender, StringListEventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewRnaSeqFastq(sender, e)));
            }
            else
            {
                foreach (var uu in rnaSeqFastqCollection)
                {
                    uu.Use = false;
                }
                foreach (var newRnaSeqFastqData in e.StringList)
                    rnaSeqFastqCollection.Add(new RNASeqFastqDataGrid(newRnaSeqFastqData));
                UpdateOutputFolderTextbox();
            }
        }

        private void MenuItem_Setup_Click(object sender, RoutedEventArgs e)
        {
            //ManageToolsFlow.Install(Path.GetDirectoryName(Assembly.GetEntryAssembly().Location));
            Spritz.Main(new string[] { "CMD.exe", "-c", "setup" });
            return;
        }

        private void MenuItem_DataDownload_Click(object sender, RoutedEventArgs e)
        {
            Spritz.Main(new string[] { "CMD.exe", "-c", "setup" });
        }

        #endregion Private Methods - Controlers

        #region Private Methods - no Controlers

        private void UpdateTaskGuiStuff()
        {
            if (staticTasksObservableCollection.Count == 0)
            {
                RunWorkflowButton.IsEnabled = false;
                RemoveLastTaskButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                ResetTasksButton.IsEnabled = false;
            }
            else
            {
                RunWorkflowButton.IsEnabled = true;
                RemoveLastTaskButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
            }
        }

        private void AddTaskToCollection(Options ye)
        {
            PreRunTask te = new PreRunTask(ye);
            staticTasksObservableCollection.Add(te);
            staticTasksObservableCollection.Last().DisplayName = "Task" + (staticTasksObservableCollection.IndexOf(te) + 1);
        }

        private void UpdateOutputFolderTextbox()
        {
            if (rnaSeqFastqCollection.Any())
            {
                var MatchingChars =
                    from len in Enumerable.Range(0, rnaSeqFastqCollection.Select(b => b.FilePath).Min(s => s.Length)).Reverse()
                    let possibleMatch = rnaSeqFastqCollection.Select(b => b.FilePath).First().Substring(0, len)
                    where rnaSeqFastqCollection.Select(b => b.FilePath).All(f => f.StartsWith(possibleMatch, StringComparison.Ordinal))
                    select possibleMatch;

                OutputFolderTextBox.Text = Path.Combine(Path.GetDirectoryName(MatchingChars.First()));
            }
            else
            {
                OutputFolderTextBox.Clear();
            }
        }

        private void AddAFile(string filepath)
        {
            var theExtension = Path.GetExtension(filepath).ToLowerInvariant();
            switch (theExtension)
            {
                case ".fa":
                    GenomeFastaDataGrid genomeFasta = new GenomeFastaDataGrid(filepath);
                    genomeFastaCollection.Add(genomeFasta);
                    break;

                case ".gtf":
                case ".gff3":
                    GeneSetDataGrid geneSet = new GeneSetDataGrid(filepath);
                    geneSetCollection.Add(geneSet);
                    break;

                case ".fastq":
                case ".fastq.gz":
                    RNASeqFastqDataGrid rnaSeqFastq = new RNASeqFastqDataGrid(filepath);
                    rnaSeqFastqCollection.Add(rnaSeqFastq);
                    UpdateOutputFolderTextbox();
                    break;
            }
        }

        #endregion Private Methods - no Controlers
    }
}