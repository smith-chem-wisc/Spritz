using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.IO;
using System.Collections.ObjectModel;
using WorkflowLayer;


namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<GenomeFastaDataGrid> genomeFastaCollection = new ObservableCollection<GenomeFastaDataGrid>();
        private readonly ObservableCollection<GeneSetDataGrid> geneSetCollection = new ObservableCollection<GeneSetDataGrid>();
        private readonly ObservableCollection<RNASeqFastqDataGrid> rnaSeqFastqCollection = new ObservableCollection<RNASeqFastqDataGrid>();
        private ObservableCollection<InRunTask> dynamicTasksObservableCollection;
        private readonly ObservableCollection<PreRunTask> staticTasksObservableCollection = new ObservableCollection<PreRunTask>();


        public MainWindow()
        {
            InitializeComponent();

            dataGridFASTA.DataContext = genomeFastaCollection;
            dataGridGeneSet.DataContext = geneSetCollection;
            dataGridRnaSeqFastq.DataContext = rnaSeqFastqCollection;
            workflowTreeView.DataContext = staticTasksObservableCollection;

            EverythingRunnerEngine.NewRnaSeqFastqHandler += AddNewRnaSeqFastq;
        }

        

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
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            
        }

        private void MenuItem_Help_Click(object sender, RoutedEventArgs e)
        {

        }

        private void MenuItem_Contact_Click(object sender, RoutedEventArgs e)
        {

        }

        private void BtnAddFastq2Proteins_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new Fastq2ProteinsWF();
            if (dialog.ShowDialog()==true)
            {

            }
        }

        private void BtnAddLncRNADiscover_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new LncRNADiscoverWFWindows();
            if (dialog.ShowDialog() == true)
            {

            }
        }

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

        private void RunWorkflowButton_Click(object sender, RoutedEventArgs e)
        {
            dynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
            for (int i = 0; i < staticTasksObservableCollection.Count; i++)
            {
                dynamicTasksObservableCollection.Add(new InRunTask("Workflow" + (i + 1) + "-" + staticTasksObservableCollection[i].spritzWorkflow.WorkflowType.ToString(), staticTasksObservableCollection[i].spritzWorkflow));
            }
            workflowTreeView.DataContext = dynamicTasksObservableCollection;
            EverythingRunnerEngine a = new EverythingRunnerEngine(dynamicTasksObservableCollection.Select(b => new Tuple<string, SpritzWorkflow>(b.DisplayName, b.workflow)).ToList(), rnaSeqFastqCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(), genomeFastaCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(), geneSetCollection.Where(b => b.Use).Select(b => b.FilePath).ToList(), OutputFolderTextBox.Text);
            var t = new Task(a.Run);
            t.Start();
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
                Filter = "Database Files|*.fa",
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
                Filter = "Database Files|*.fastq",
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
                Filter = "Database Files|*.gtf;*.gff3",
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

        private void LoadTaskButton_Click(object sender, RoutedEventArgs e)
        {

        }

        private void ClearTasksButton_Click(object sender, RoutedEventArgs e)
        {

        }

        private void RemoveLastTaskButton_Click(object sender, RoutedEventArgs e)
        {

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

        private void UpdateOutputFolderTextbox()
        {
            if (rnaSeqFastqCollection.Any())
            {
                var MatchingChars =
                    from len in Enumerable.Range(0, rnaSeqFastqCollection.Select(b => b.FilePath).Min(s => s.Length)).Reverse()
                    let possibleMatch = rnaSeqFastqCollection.Select(b => b.FilePath).First().Substring(0, len)
                    where rnaSeqFastqCollection.Select(b => b.FilePath).All(f => f.StartsWith(possibleMatch, StringComparison.Ordinal))
                    select possibleMatch;

                OutputFolderTextBox.Text = Path.Combine(Path.GetDirectoryName(MatchingChars.First()), @"$DATETIME");
            }
            else
            {
                OutputFolderTextBox.Clear();
            }
        }
    }
}
