﻿using System;
using System.Collections.ObjectModel;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Threading;
using SpritzBackend;

namespace Spritz
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private readonly ObservableCollection<RNASeqFastqDataGrid> RnaSeqFastqCollection = new();
        private ObservableCollection<InRunTask> DynamicTasksObservableCollection = new();
        private readonly ObservableCollection<PreRunTask> StaticTasksObservableCollection = new();
        private readonly ObservableCollection<SRADataGrid> SraCollection = new();
        private RunnerEngine Runner;
        private readonly Regex outputScrub = new(@"(\[\d+m)");

        public int DockerCPUs { get; set; }
        public double DockerMemory { get; set; }
        private string DockerImage { get; set; } = "smithlab/spritz";
        private string DockerStdOut { get; set; }
        private bool ShowStdOut { get; set; } = true;
        private string DockerSystemInfo { get; set; }
        private bool IsRunning { get; set; }

        public MainWindow()
        {
            try
            {
                InitializeComponent();
                DataGridRnaSeqFastq.DataContext = RnaSeqFastqCollection;
                WorkflowTreeView.DataContext = StaticTasksObservableCollection;
                LbxSRAs.ItemsSource = SraCollection;

                // Version information
                try
                {
                    SpritzUpdater.GetVersionNumbersFromWeb();
                }
                catch (Exception e)
                {
                    MessageBox.Show("Could not get newest version from web: " + e.Message, "Setup", MessageBoxButton.OK, MessageBoxImage.Warning);
                }

                // Check Docker setup
                Dispatcher.Invoke(() =>
                {
                    Process proc = new();
                    proc.StartInfo.FileName = "Powershell.exe";
                    proc.StartInfo.Arguments = "docker system info";
                    proc.StartInfo.UseShellExecute = false;
                    proc.StartInfo.RedirectStandardOutput = true;
                    proc.StartInfo.RedirectStandardError = true;
                    proc.StartInfo.CreateNoWindow = true;
                    proc.Start();
                    StreamReader outputReader = proc.StandardOutput;
                    DockerSystemInfo = outputReader.ReadToEnd();
                    proc.WaitForExit();
                });
                bool isDockerInstalled = !string.IsNullOrEmpty(DockerSystemInfo) && !DockerSystemInfo.Contains("error", StringComparison.InvariantCultureIgnoreCase);
                if (isDockerInstalled)
                {
                    ParseDockerSystemInfo(DockerSystemInfo);
                }
                string message = isDockerInstalled ?
                    "In Docker Desktop, please ensure all shared drives are enabled, and please ensure a Disk image size of at least 80 GB is enabled." :
                    "Docker is not running or not installed. Please have Docker Desktop installed and running, enable all shared drives, and ensure a Disk image size of at least 80 GB is enabled.";
                if (isDockerInstalled && DockerMemory < 16)
                {
                    message += $"{Environment.NewLine}{Environment.NewLine}The memory allocated to Docker is low ({DockerMemory}GB). Please raise this value above 16 GB in Docker Desktop if possible.";
                }
                MessageBox.Show(message, "Setup", MessageBoxButton.OK, isDockerInstalled ? MessageBoxImage.Information : MessageBoxImage.Error);
            }
            catch (Exception e)
            {
                MessageBox.Show($"ERROR: {e.Message}", "Initialization Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void ParseDockerSystemInfo(string dockerSystemInfo)
        {
            string[] infoLines = dockerSystemInfo.Split('\n');
            string cpuLine = infoLines.FirstOrDefault(line => line.Trim().StartsWith("CPUs", StringComparison.InvariantCultureIgnoreCase));
            if (cpuLine != null && int.TryParse(cpuLine.Split(':')[1].Trim(), out int dockerThreads))
            {
                DockerCPUs = dockerThreads;
            }

            double gibToGbConversion = 1.07374;
            string memoryLine = infoLines.FirstOrDefault(line => line.Trim().StartsWith("Total Memory", StringComparison.InvariantCultureIgnoreCase));
            if (memoryLine != null && double.TryParse(memoryLine.Split(':')[1].Replace("GiB", "").Trim(), out double memoryGB))
            {
                DockerMemory = memoryGB * gibToGbConversion;
            }
        }

        protected override void OnClosed(EventArgs e)
        {
            string message = "Are you sure you would like to exit Spritz?";
            message += IsRunning ? " This will stop all Spritz processes, which may take a few moments." : "";
            if (MessageBox.Show(message, "Exit Spritz", MessageBoxButton.OKCancel) == MessageBoxResult.Cancel)
                return;
            StopDocker("stop"); // may need to kill processes if we see that they get stuck in the future, but killing leaves some files open, which makes them hard to delete
            base.OnClosed(e);
        }

        private void CancelTasksButton_Click(object sender, RoutedEventArgs e)
        {
            StopDocker("stop");
        }

        private void StopDocker(string command)
        {
            // new process that kills docker container (if any)
            if (Runner != null && !string.IsNullOrEmpty(Runner.PathToWorkflow))
            {
                Process proc = new();
                proc.StartInfo.FileName = "Powershell.exe";
                proc.StartInfo.Arguments = $"docker {command} {Runner.SpritzContainerName}";
                proc.StartInfo.CreateNoWindow = true;
                proc.StartInfo.UseShellExecute = false;
                proc.Start();

                if (proc != null && !proc.HasExited)
                {
                    proc.WaitForExit();
                }
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
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            if (SpritzUpdater.NewestKnownVersion != null && !SpritzUpdater.IsVersionLower(SpritzUpdater.NewestKnownVersion))
            {
                try
                {
                    SpritzUpdater newwind = new();
                    newwind.ShowDialog();
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.ToString());
                }
            }
        }

        private void MenuItem_Wiki_Click(object sender, RoutedEventArgs e)
        {
            Process.Start(@"https://github.com/smith-chem-wisc/Spritz/wiki");
        }

        private void MenuItem_Contact_Click(object sender, RoutedEventArgs e)
        {
            Process.Start(@"https://github.com/smith-chem-wisc/Spritz");
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

                DynamicTasksObservableCollection = new();
                DynamicTasksObservableCollection.Add(new InRunTask("Workflow 1", StaticTasksObservableCollection.First().options));
                WorkflowTreeView.DataContext = DynamicTasksObservableCollection;

                Runner = new RunnerEngine(DynamicTasksObservableCollection.Select(b => new Tuple<string, SpritzOptions>(b.DisplayName, b.options)).First(), OutputFolderTextBox.Text);

                InformationTextBox.Document.Blocks.Clear();
                Runner.GenerateSpritzCMDCommand(StaticTasksObservableCollection.First().options);
                InformationTextBox.AppendText($"Command executing: Powershell.exe {Runner.GenerateCommandsDry(DockerImage, Runner.SpritzCMDCommand)}\n\n"); // keep for debugging
                InformationTextBox.AppendText($"Saving output to {Runner.PathToWorkflow}. Please monitor it there...\n\n");

                IsRunning = true;
                var t = new Task(RunEverythingRunner);
                t.Start();
                t.ContinueWith(DisplayAnyErrors);

                // update gui
                RunWorkflowButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = true;
                BtnWorkFlow.IsEnabled = false;
                ResetTasksButton.IsEnabled = true;
            }
            catch (TaskCanceledException)
            {
                // Ignore error
            }
        }

        private void RunEverythingRunner()
        {
            Process proc = new();
            proc.StartInfo.FileName = "Powershell.exe";
            proc.StartInfo.Arguments = Runner.GenerateCommandsDry(DockerImage, Runner.SpritzCMDCommand);
            proc.StartInfo.UseShellExecute = false;
            proc.StartInfo.RedirectStandardOutput = true;
            proc.StartInfo.RedirectStandardError = true;
            proc.StartInfo.CreateNoWindow = true;
            proc.OutputDataReceived += new DataReceivedEventHandler(OutputHandler);
            proc.ErrorDataReceived += new DataReceivedEventHandler(OutputHandler);
            proc.Start();
            proc.BeginOutputReadLine();
            proc.BeginErrorReadLine();
            proc.WaitForExit();
        }

        private void OutputHandler(object source, DataReceivedEventArgs e)
        {
            Dispatcher.Invoke(() =>
            {
                if (!string.IsNullOrEmpty(e.Data))
                {
                    string output = outputScrub.Replace(e.Data, "");
                    DockerStdOut += output + Environment.NewLine;
                    if (ShowStdOut)
                    {
                        lock (InformationTextBox)
                        {
                            InformationTextBox.AppendText(output + Environment.NewLine);
                        }
                    }
                    using StreamWriter sw = File.Exists(Runner.PathToWorkflow) ? File.AppendText(Runner.PathToWorkflow) : File.CreateText(Runner.PathToWorkflow);
                    sw.WriteLine(output);
                }
            });
        }

        private void DisplayAnyErrors(Task obj)
        {
            Dispatcher.Invoke(() => InformationTextBox.AppendText("Done!" + Environment.NewLine));
            if (StaticTasksObservableCollection.Count > 0)
            {
                Dispatcher.Invoke(() => MessageBox.Show("Finished! Workflow summary is located in "
                    + StaticTasksObservableCollection.First().options.AnalysisDirectory, "Spritz Workflow",
                    MessageBoxButton.OK, MessageBoxImage.Information));
            }
            IsRunning = false;
        }

        private void BtnAddRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new()
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
        }

        private void BtnClearRnaSeqFastq_Click(object sender, RoutedEventArgs e)
        {
            RnaSeqFastqCollection.Clear();
            UpdateOutputFolderTextbox();
        }

        private void ClearTasksButton_Click(object sender, RoutedEventArgs e)
        {
            StaticTasksObservableCollection.Clear();
            WorkflowTreeView.DataContext = StaticTasksObservableCollection;
            InformationTextBox.Document.Blocks.Clear();
            UpdateTaskGuiStuff();
        }

        private void ResetTasksButton_Click(object sender, RoutedEventArgs e)
        {
            RunWorkflowButton.IsEnabled = true;
            ClearTasksButton.IsEnabled = true;
            BtnWorkFlow.IsEnabled = false;
            ResetTasksButton.IsEnabled = false;

            DynamicTasksObservableCollection.Clear();
            WorkflowTreeView.DataContext = StaticTasksObservableCollection;
        }

        private void BtnAddSRA_Click(object sender, RoutedEventArgs e)
        {
            AddSRA(TbxSRA.Text, true);
        }

        private void BtnAddSRASingleEnd_Click(object sender, RoutedEventArgs e)
        {
            AddSRA(TbxSRA.Text, false);
        }

        private void AddSRA(string name, bool isPairedEnd)
        {
            if (name.Contains("SR") || name.Contains("ER"))
            {
                if (SraCollection.Any(s => s.Name == name.Trim()))
                {
                    MessageBox.Show("That SRA has already been added. Please choose a new SRA accession.", "Workflow", MessageBoxButton.OK, MessageBoxImage.Information);
                }
                else
                {
                    SRADataGrid sraDataGrid = new(name.Trim(), isPairedEnd);
                    SraCollection.Add(sraDataGrid);
                }
            }
            else if (MessageBox.Show("SRA accessions are expected to start with \"SR\" or \"ER\", such as SRX254398 or ERR315327. View the GEO SRA website?", "Workflow", MessageBoxButton.YesNo, MessageBoxImage.Question, MessageBoxResult.No) == MessageBoxResult.Yes)
            {
                Process.Start("https://www.ncbi.nlm.nih.gov/sra");
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
                var dialog = new WorkFlowWindow(string.IsNullOrEmpty(OutputFolderTextBox.Text) ? 
                    SpritzOptions.DefaultAnalysisDirectory() : 
                    OutputFolderTextBox.Text);

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

        private void UpdateTaskGuiStuff()
        {
            if (StaticTasksObservableCollection.Count == 0)
            {
                RunWorkflowButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
                BtnWorkFlow.IsEnabled = true;
                ResetTasksButton.IsEnabled = false;
            }
            else
            {
                RunWorkflowButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
                BtnWorkFlow.IsEnabled = false;
                ResetTasksButton.IsEnabled = false;
            }
        }

        private void AddTaskToCollection(SpritzOptions ye)
        {
            PreRunTask te = new(ye);
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
            var theExtension = Path.GetExtension(filepath).ToLowerInvariant();
            string basepath = theExtension == ".gz" ? Path.GetFileNameWithoutExtension(Path.GetFileNameWithoutExtension(filepath)) : Path.GetFileNameWithoutExtension(filepath);
            theExtension = theExtension == ".gz" ? Path.GetExtension(Path.GetFileNameWithoutExtension(filepath)).ToLowerInvariant() : theExtension;
            switch (theExtension)
            {
                case ".fastq":
                    if (Path.GetFileName(basepath).EndsWith("_1") || Path.GetFileName(basepath).EndsWith("_2"))
                    {
                        bool isPairedEnd = Path.GetFileName(basepath).EndsWith("_2");
                        bool isUnpairedSecondMate = false;
                        bool isDifferentFolder = false;
                        foreach (var fq in RnaSeqFastqCollection)
                        {
                            bool hasMatePair = Path.GetFileName(fq.FileName.Substring(0, fq.FileName.Length - 2)) ==
                                Path.GetFileName(basepath.Substring(0, basepath.Length - 2));
                            if (hasMatePair)
                                fq.IsPairedEnd = true;
                            isPairedEnd = isPairedEnd || hasMatePair;
                            isUnpairedSecondMate = isUnpairedSecondMate || hasMatePair;
                            isDifferentFolder = isDifferentFolder ||
                                Path.GetDirectoryName(fq.FilePath) != Path.GetDirectoryName(filepath);
                        }
                        if (isDifferentFolder)
                        {
                            MessageBox.Show("All FASTQs must be located in the same analysis directory.", "Add FASTQ", MessageBoxButton.OK, MessageBoxImage.Warning);
                        }
                        else
                        {
                            RNASeqFastqDataGrid rnaSeqFastq = new(filepath, isPairedEnd);
                            RnaSeqFastqCollection.Add(rnaSeqFastq);
                            UpdateOutputFolderTextbox();
                        }
                        break;
                    }
                    else
                    {
                        MessageBox.Show("FASTQ files must have *_1.fastq and/or *_2.fastq extensions.\n\n" +
                            "Single-end experiments must have unique filenames preceeding the *_1.fastq ending.\n\n" +
                            "Paired-end experiments must have the same unique filename preceeding the *_1.fastq and *_2.fastq suffixes.",
                            "Run Workflows", MessageBoxButton.OK, MessageBoxImage.Information);
                        return;
                    }
            }
        }

        private void WorkflowTreeView_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            //var a = sender as TreeView;
            //if (a.SelectedItem is PreRunTask preRunTask)
            //{
            //    var workflowDialog = new WorkFlowWindow(preRunTask.options);
            //    workflowDialog.ShowDialog();
            //    WorkflowTreeView.Items.Refresh();
            //}
        }

        private void WarningsTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            InformationTextBox.ScrollToEnd();
        }

        private void DockerImage_TextChanged(object sender, TextChangedEventArgs e)
        {
            DockerImage = tb_DockerImage.Text;
        }

        private void ShowTopButton_Click(object sender, RoutedEventArgs e)
        {
            ShowStdOut = false;
            Dispatcher.Invoke(() =>
            {
                InformationTextBox.Document.Blocks.Clear();

                Process proc = new();
                proc.StartInfo.FileName = "Powershell.exe";
                proc.StartInfo.Arguments = Runner.GenerateTopComand();
                proc.StartInfo.UseShellExecute = false;
                proc.StartInfo.RedirectStandardOutput = true;
                proc.StartInfo.CreateNoWindow = true;
                proc.Start();
                StreamReader outputReader = proc.StandardOutput;
                InformationTextBox.AppendText(outputReader.ReadToEnd());
                proc.WaitForExit();
            });
        }

        private void ShowOutputButton_Click(object sender, RoutedEventArgs e)
        {
            ShowStdOut = true;
            lock (InformationTextBox)
            {
                InformationTextBox.Document.Blocks.Clear();
                InformationTextBox.AppendText(DockerStdOut);
            }
        }
    }
}
