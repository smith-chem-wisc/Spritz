using System;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RNASeqFastqDataGrid> RnaSeqFastqCollection = new ObservableCollection<RNASeqFastqDataGrid>();
        private ObservableCollection<InRunTask> DynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
        private readonly ObservableCollection<PreRunTask> StaticTasksObservableCollection = new ObservableCollection<PreRunTask>();
        private readonly ObservableCollection<SRADataGrid> SraCollection = new ObservableCollection<SRADataGrid>();
        private CancellationTokenSource TokenSource = new CancellationTokenSource();
        private EverythingRunnerEngine Everything;
        //private Task EverythingTask;

        public MainWindow()
        {
            InitializeComponent();
            DataGridRnaSeqFastq.DataContext = RnaSeqFastqCollection;
            workflowTreeView.DataContext = StaticTasksObservableCollection;
            LbxSRAs.ItemsSource = SraCollection;
            MessageBox.Show("Please have Docker Desktop installed. Under \"Shared Drives\", select drives to be shared and click \"Apply\".", "Setup", MessageBoxButton.OK, MessageBoxImage.Information);
        }

        protected override void OnClosed(EventArgs e)
        {
            // TODO: implement some way of killing EverythingTask

            // new process that kills docker container (if any)
            Process proc = new Process();
            proc.StartInfo.FileName = "Powershell.exe";
            proc.StartInfo.Arguments = "docker kill spritz";
            proc.StartInfo.CreateNoWindow = true;
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardError = true;
            proc.Start();

            if (proc != null && !proc.HasExited)
            {
                proc.WaitForExit();
            }

            base.OnClosed(e);
        }

        private void UpdateSRABox()
        {
            if (RnaSeqFastqCollection.Count > 0)
            {
                TbxSRA.IsEnabled = false;
                BtnAddSRA.IsEnabled = false;
                BtnClearSRA.IsEnabled = false;
            }
            else
            {
                TbxSRA.IsEnabled = true;
                BtnAddSRA.IsEnabled = true;
                BtnClearSRA.IsEnabled = true;
            }
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
                    DataGridRnaSeqFastq.Items.Refresh();
                }
            }
            UpdateOutputFolderTextbox();
            UpdateSRABox();
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

        private void RunWorkflowButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                if (StaticTasksObservableCollection.Count == 0)
                {
                    MessageBox.Show("You must add a workflow before a run.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }
                else if (RnaSeqFastqCollection.Any() && GetPathToFastqs().CompareTo(OutputFolderTextBox.Text) != 0) // to be edited
                {
                    MessageBox.Show("FASTQ files do not exist in the user-defined analysis directory.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                    return;
                }

                DynamicTasksObservableCollection = new ObservableCollection<InRunTask>();
                for (int i = 0; i < StaticTasksObservableCollection.Count; i++)
                {
                    DynamicTasksObservableCollection.Add(new InRunTask("Workflow" + (i + 1), StaticTasksObservableCollection[i].options));
                }
                workflowTreeView.DataContext = DynamicTasksObservableCollection;
                Everything = new EverythingRunnerEngine(DynamicTasksObservableCollection.Select(b => new Tuple<string, Options>(b.DisplayName, b.options)).ToList(), OutputFolderTextBox.Text);
                //WarningsTextBox.AppendText(string.Join("\n", Everything.GenerateCommandsDry().Select(x => $"Command executing: CMD.exe {x}"))); // keep for debugging
                WarningsTextBox.AppendText(string.Join("\n", Everything.GenerateCommandsDry().Select(x => $"Command executing: Powershell.exe {x}"))); // keep for debugging
                var t = new Task(Everything.Run);
                t.Start();
                t.ContinueWith(DisplayAnyErrors);
                RunWorkflowButton.IsEnabled = false;
            }
            catch (TaskCanceledException)
            {
                // Ignore error
            }
        }

        private void DisplayAnyErrors(Task obj)
        {
            if (Everything.StdErr != null && Everything.StdErr != "")
            {
                var message = "Run failed, Exception: " + Everything.StdErr;
                Dispatcher.Invoke(() => WarningsTextBox.AppendText(message + Environment.NewLine));
                var messageBoxResult = MessageBox.Show(message + "\n\nWould you like to report this crash?", "Runtime Error", MessageBoxButton.YesNo, MessageBoxImage.Question);
                if (messageBoxResult == MessageBoxResult.Yes)
                {
                    string body = Everything.StdErr;
                    //+ "%0D%0A" + exception.Data +
                    //"%0D%0A" + exception.StackTrace +
                    //"%0D%0A" + exception.Source +
                    //"%0D%0A %0D%0A %0D%0A %0D%0A SYSTEM INFO: %0D%0A ???" + // TODO: implement this system info check
                    //"%0D%0A%0D%0A Spritz: version ???" + // TODO: implement this version check.
                    //"%0D%0A %0D%0A %0D%0A %0D%0A TOML: %0D%0A " +
                    //tomlText;
                    body = body.Replace('&', ' ');
                    body = body.Replace("\n", "%0D%0A");
                    body = body.Replace("\r", "%0D%0A");
                    string mailto = string.Format("mailto:{0}?Subject=Spritz. Issue:&Body={1}", "mm_support@chem.wisc.edu", body);
                    System.Diagnostics.Process.Start(mailto);
                    Console.WriteLine(body);
                }
            }
            else
            {
                Dispatcher.Invoke(() => WarningsTextBox.AppendText("Done!" + Environment.NewLine));
                Dispatcher.Invoke(() => MessageBox.Show("Finished!", "Spritz Workflow", MessageBoxButton.OK, MessageBoxImage.Information));
            }
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
            {
                foreach (var filepath in openPicker.FileNames)
                {
                    AddAFile(filepath);
                }
            }
            DataGridRnaSeqFastq.Items.Refresh();
            UpdateSRABox();
        }

        private void BtnClearRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            RnaSeqFastqCollection.Clear();
            UpdateOutputFolderTextbox();
            UpdateSRABox();
        }

        //private void LoadTaskButton_Click(object sender, RoutedEventArgs e)
        //{
        //    Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
        //    {
        //        Filter = "TOML files(*.toml)|*.toml",
        //        FilterIndex = 1,
        //        RestoreDirectory = true,
        //        Multiselect = true
        //    };
        //    if (openPicker.ShowDialog() == true)
        //    {
        //        foreach (var tomlFromSelected in openPicker.FileNames)
        //        {
        //            AddAFile(tomlFromSelected);
        //        }
        //    }
        //    UpdateTaskGuiStuff();
        //}

        private void ClearTasksButton_Click(object sender, RoutedEventArgs e)
        {
            StaticTasksObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        //private void RemoveLastTaskButton_Click(object sender, RoutedEventArgs e)
        //{
        //    StaticTasksObservableCollection.RemoveAt(StaticTasksObservableCollection.Count - 1);
        //    UpdateTaskGuiStuff();
        //}

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
            RunWorkflowButton.IsEnabled = true;
            ResetTasksButton.IsEnabled = false;
            for (int i = 0; i < DynamicTasksObservableCollection.Count; i++)
            {
                StaticTasksObservableCollection.Add(new PreRunTask(StaticTasksObservableCollection[i].options));
            }
            DynamicTasksObservableCollection.Clear();
            workflowTreeView.DataContext = StaticTasksObservableCollection;
        }

        //private void AddNewRnaSeqFastq(object sender, StringListEventArgs e)
        //{
        //    if (!Dispatcher.CheckAccess())
        //    {
        //        Dispatcher.BeginInvoke(new Action(() => AddNewRnaSeqFastq(sender, e)));
        //    }
        //    else
        //    {
        //        foreach (var uu in RnaSeqFastqCollection)
        //        {
        //            uu.Use = false;
        //        }
        //        foreach (var newRnaSeqFastqData in e.StringList)
        //        {
        //            RnaSeqFastqCollection.Add(new RNASeqFastqDataGrid(newRnaSeqFastqData));
        //        }
        //        UpdateOutputFolderTextbox();
        //    }
        //}

        private void BtnAddSRA_Click(object sender, RoutedEventArgs e)
        {
            if (TbxSRA.Text.Contains("SR") || TbxSRA.Text.Contains("ER"))
            {
                if (SraCollection.Any(s => s.Name == TbxSRA.Text))
                {
                    MessageBox.Show("That SRA has already been added. Please choose a new SRA accession.", "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
                }
                else
                {
                    SRADataGrid sraDataGrid = new SRADataGrid(TbxSRA.Text.Trim());
                    SraCollection.Add(sraDataGrid);
                }
            }
            else if (MessageBox.Show("SRA accessions are expected to start with \"SR\" or \"ER\", such as SRX254398 or ERR315327. View the GEO SRA website?", "Workflow", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.No) == MessageBoxResult.Yes)
            {
                System.Diagnostics.Process.Start("https://www.ncbi.nlm.nih.gov/sra");
            }
        }

        private void BtnClearSRA_Click(object sender, RoutedEventArgs e)
        {
            SraCollection.Clear();
            BtnAddSRA.IsEnabled = true;
        }

        private void BtnWorkFlow_Click(object sender, RoutedEventArgs e)
        {
            if (SraCollection.Count == 0 && RnaSeqFastqCollection.Count == 0)
            {
                if (MessageBox.Show("You have not added any nucleic acid sequencing data (SRA accession or fastq files). Would you like to continue to make a protein database from the reference gene model?", "Workflow", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.No) == MessageBoxResult.No)
                {
                    return;
                }
            }

            try
            {
                var dialog = new WorkFlowWindow(OutputFolderTextBox.Text == "" ? new Options().AnalysisDirectory : OutputFolderTextBox.Text);
                if (dialog.ShowDialog() == true)
                {
                    AddTaskToCollection(dialog.Options);
                    UpdateTaskGuiStuff();
                    UpdateOutputFolderTextbox();
                }
            }
            catch (InvalidOperationException)
            {
                // does not open workflow window until all fastq files are added, if any
            }
        }

        //private void BtnSaveRnaSeqFastqSet_Click(object sender, RoutedEventArgs e)
        //{
        //    try
        //    {
        //        WriteExperDesignToTsv(OutputFolderTextBox.Text);
        //    }
        //    catch (Exception ex)
        //    {
        //        MessageBox.Show("Could not save experimental design!\n\n" + ex.Message, "Experimental Design", MessageBoxButton.OK, MessageBoxImage.Warning);
        //        return;
        //    }
        //}

        private void UpdateTaskGuiStuff()
        {
            if (StaticTasksObservableCollection.Count == 0)
            {
                RunWorkflowButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                BtnWorkFlow.IsEnabled = true;
            }
            else
            {
                RunWorkflowButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
                BtnWorkFlow.IsEnabled = false;
            }
        }

        private void AddTaskToCollection(Options ye)
        {
            PreRunTask te = new PreRunTask(ye);
            StaticTasksObservableCollection.Add(te);
            StaticTasksObservableCollection.Last().DisplayName = "Task" + (StaticTasksObservableCollection.IndexOf(te) + 1);
        }

        private string GetPathToFastqs()
        {
            var MatchingChars =
                    from len in Enumerable.Range(0, RnaSeqFastqCollection.Select(b => b.FilePath).Min(s => s.Length)).Reverse()
                    let possibleMatch = RnaSeqFastqCollection.Select(b => b.FilePath).First().Substring(0, len)
                    where RnaSeqFastqCollection.Select(b => b.FilePath).All(f => f.StartsWith(possibleMatch, StringComparison.Ordinal))
                    select possibleMatch;

            return Path.Combine(Path.GetDirectoryName(MatchingChars.First()));
        }

        //private string GetPathToFastqDirectory(string path)
        //{
        //    var filePath = path.Split('\\');
        //    var newPath = string.Join("\\", filePath.Take(filePath.Length - 1));
        //    return newPath;
        //}

        //private void UpdateOutputFolderTextbox(string filePath = null)
        //{
        //    // if new files have different path than current text in output, then throw error
        //    if (StaticTasksObservableCollection.Count > 0)
        //    {
        //        OutputFolderTextBox.Text = StaticTasksObservableCollection.First().options.AnalysisDirectory;
        //    }
        //    else if (RnaSeqFastqCollection.Any())
        //    {
        //        if (filePath != null && OutputFolderTextBox.Text != "" && GetPathToFastqDirectory(filePath).CompareTo(OutputFolderTextBox.Text) != 0)
        //        {
        //            throw new InvalidOperationException();
        //        }
        //        OutputFolderTextBox.Text = GetPathToFastqs();
        //    }
        //    else
        //    {
        //        OutputFolderTextBox.Clear();
        //    }
        //}

        private void UpdateOutputFolderTextbox()
        {
            if (StaticTasksObservableCollection.Count > 0)
            {
                OutputFolderTextBox.Text = StaticTasksObservableCollection.First().options.AnalysisDirectory;
            }
            else if (RnaSeqFastqCollection.Any())
            {
                OutputFolderTextBox.Text = GetPathToFastqs();
            }
            else
            {
                OutputFolderTextBox.Clear();
            }
        }

        private void AddAFile(string filepath)
        {
            if (SraCollection.Count == 0)
            {
                var theExtension = Path.GetExtension(filepath).ToLowerInvariant();
                theExtension = theExtension == ".gz" ? Path.GetExtension(Path.GetFileNameWithoutExtension(filepath)).ToLowerInvariant() : theExtension;
                switch (theExtension)
                {
                    case ".fastq":
                        if (Path.GetFileName(filepath).Contains("_1") || Path.GetFileName(filepath).Contains("_2"))
                        {
                            RNASeqFastqDataGrid rnaSeqFastq = new RNASeqFastqDataGrid(filepath);
                            RnaSeqFastqCollection.Add(rnaSeqFastq);
                            UpdateOutputFolderTextbox();
                            break;
                        }
                        else
                        {
                            MessageBox.Show("FASTQ files must have *_1.fastq and *_2.fastq extensions.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                            return;
                        }

                        //case ".toml":
                        //    TomlTable tomlFile = null;
                        //    try
                        //    {
                        //        tomlFile = Toml.ReadFile(filepath);
                        //    }
                        //    catch (Exception)
                        //    {
                        //        break;
                        //    }
                        //    var ye1 = Toml.ReadFile<Options>(filepath);
                        //    AddTaskToCollection(ye1);
                        //    break;
                }
            }
            else
            {
                MessageBox.Show("User already added SRA number. Please only choose one input: 1) SRA accession 2) FASTQ files.", "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                return;
            }
        }

        //private bool InstallationDialogAndCheck()
        //{
        //    var exePath = Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location);

        //    if (!WrapperUtility.CheckToolSetup(exePath))
        //    {
        //        try
        //        {
        //            var dialog = new InstallWindow();
        //            dialog.ShowDialog();
        //        }
        //        catch (Exception ex)
        //        {
        //            MessageBox.Show(ex.ToString(), "Installation", MessageBoxButton.OK, MessageBoxImage.Error);
        //        }
        //    }
        //    return WrapperUtility.CheckBashSetup() && WrapperUtility.CheckToolSetup(Environment.CurrentDirectory);
        //}

        private void WriteExperDesignToTsv(string filePath)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                foreach (var aFastq in RnaSeqFastqCollection)
                {
                    output.WriteLine(aFastq.FileName +
                        "\t" + aFastq.Experiment +
                        "\t" + aFastq.MatePair);
                }
            }
        }

        private void workflowTreeView_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var a = sender as TreeView;
            if (a.SelectedItem is PreRunTask preRunTask)
            {
                var workflowDialog = new WorkFlowWindow(preRunTask.options);
                workflowDialog.ShowDialog();
                workflowTreeView.Items.Refresh();
            }
        }
    }
}