using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using Spritz.Avalonia.Services;
using SpritzBackend;

namespace Spritz.Avalonia.ViewModels;

/// <summary>
/// Everything the window does, with no UI types in it.
///
/// The WPF version keeps this logic in MainWindow.xaml.cs, so none of it can be tested without a
/// window. Here the view binds to these properties and commands, which means the run flow, the
/// docker detection and the argument construction are all reachable from a test.
/// </summary>
internal sealed partial class MainWindowViewModel : ObservableObject
{
    private static readonly Regex AnsiEscape = new(@"(\[\d+m)", RegexOptions.Compiled);

    private readonly DockerProcessRunner _docker;
    private CancellationTokenSource _cancellation;
    private RunnerEngine _runner;

    public ObservableCollection<SRADataGrid> Sras { get; } = new();
    public ObservableCollection<RNASeqFastqDataGrid> Fastqs { get; } = new();
    public ObservableCollection<PreRunTask> Workflows { get; } = new();

    [ObservableProperty] private string _sraAccession = string.Empty;
    [ObservableProperty] private string _dockerImage = "smithlab/spritz";
    [ObservableProperty] private string _outputFolder = string.Empty;
    [ObservableProperty] private string _information = string.Empty;
    [ObservableProperty] private string _dockerStatus = "Checking for Docker...";
    [ObservableProperty] private bool _isRunning;

    public MainWindowViewModel()
    {
        _docker = new DockerProcessRunner(AppendInformation);
        _ = CheckDockerAsync();
    }

    /// <summary>
    /// Reports whether Docker is usable, and how much memory it has. The WPF version does this in
    /// the constructor with a blocking Dispatcher.Invoke and then a modal dialog; here it is
    /// asynchronous and surfaces as a status line, so a slow or absent Docker cannot hang startup.
    /// </summary>
    private async Task CheckDockerAsync()
    {
        var (available, info) = await _docker.QuerySystemInfoAsync();
        if (!available)
        {
            DockerStatus = "Docker is not running or not installed. Install Docker and start it, "
                + "then restart Spritz.";
            return;
        }

        double memoryGb = ParseTotalMemoryGb(info);
        DockerStatus = memoryGb > 0 && memoryGb < 16
            ? $"Docker is running, but only {memoryGb:F1} GB is allocated to it. Raise this above 16 GB if you can."
            : "Docker is running.";
    }

    /// <summary>Reads "Total Memory: 31.2GiB" out of `docker system info`, in GB.</summary>
    internal static double ParseTotalMemoryGb(string dockerSystemInfo)
    {
        const double gibToGb = 1.07374;
        string line = dockerSystemInfo
            .Split('\n')
            .FirstOrDefault(l => l.TrimStart().StartsWith("Total Memory", StringComparison.OrdinalIgnoreCase));
        if (line is null || !line.Contains(':'))
        {
            return 0;
        }
        return double.TryParse(line.Split(':')[1].Replace("GiB", string.Empty).Trim(), out double gib)
            ? gib * gibToGb
            : 0;
    }

    [RelayCommand]
    private void AddPairedEndSra() => AddSra(isPairedEnd: true);

    [RelayCommand]
    private void AddSingleEndSra() => AddSra(isPairedEnd: false);

    private void AddSra(bool isPairedEnd)
    {
        string accession = SraAccession?.Trim();
        if (string.IsNullOrEmpty(accession))
        {
            return;
        }
        if (!accession.StartsWith("SR", StringComparison.OrdinalIgnoreCase)
            && !accession.StartsWith("ER", StringComparison.OrdinalIgnoreCase))
        {
            AppendInformation($"'{accession}' does not look like an SRA accession; expected something like SRR13737862.");
            return;
        }
        if (Sras.Any(s => s.Name == accession))
        {
            return;
        }
        Sras.Add(new SRADataGrid(accession, isPairedEnd));
        SraAccession = string.Empty;
    }

    [RelayCommand]
    private void ClearSras() => Sras.Clear();

    [RelayCommand]
    private void ClearFastqs() => Fastqs.Clear();

    /// <summary>Adds FASTQ paths chosen by the view, which owns the file dialog.</summary>
    public void AddFastqFiles(IEnumerable<string> paths)
    {
        foreach (string path in paths.Where(p => p.EndsWith(".fastq") || p.EndsWith(".fastq.gz")))
        {
            if (Fastqs.All(f => f.FilePath != path))
            {
                Fastqs.Add(new RNASeqFastqDataGrid(path, isPairedEnd: path.Contains("_1") || path.Contains("_2")));
            }
        }
        if (string.IsNullOrEmpty(OutputFolder) && Fastqs.Count > 0)
        {
            OutputFolder = Path.GetDirectoryName(Fastqs[0].FilePath) ?? string.Empty;
        }
    }

    public bool CanRun => !IsRunning && Workflows.Count > 0 && !string.IsNullOrWhiteSpace(OutputFolder);

    [RelayCommand]
    private async Task RunWorkflowAsync()
    {
        if (Workflows.Count == 0)
        {
            AppendInformation("Add a workflow before running.");
            return;
        }
        if (string.IsNullOrWhiteSpace(OutputFolder))
        {
            AppendInformation("Choose an analysis directory before running.");
            return;
        }

        SpritzOptions options = Workflows.First().options;
        _runner = new RunnerEngine(new Tuple<string, SpritzOptions>("Workflow 1", options), OutputFolder);
        _cancellation = new CancellationTokenSource();
        IsRunning = true;

        try
        {
            IReadOnlyList<string> containerCommand = DockerArguments.SpritzCmd(options);
            IReadOnlyList<string> runArguments = DockerArguments.Run(
                DockerImage, _runner.SpritzContainerName, _runner.AnalysisDirectory,
                _runner.ResourcesDirectory, containerCommand);

            if (DockerArguments.ShouldPull(DockerImage))
            {
                AppendInformation(DockerArguments.Describe(DockerArguments.Pull(DockerImage)));
                await _docker.RunAsync(DockerArguments.Pull(DockerImage), _cancellation.Token);
            }

            AppendInformation(DockerArguments.Describe(runArguments));
            AppendInformation($"Saving output to {_runner.PathToWorkflow}");
            AppendInformation(string.Empty);

            int exitCode = await _docker.RunAsync(runArguments, _cancellation.Token);
            AppendInformation(exitCode == 0
                ? $"Done. Results are in {options.AnalysisDirectory}"
                : $"Spritz exited with code {exitCode}. See the output above.");
        }
        catch (Exception exception)
        {
            AppendInformation($"ERROR: {exception.Message}");
        }
        finally
        {
            IsRunning = false;
        }
    }

    [RelayCommand]
    private async Task CancelAsync()
    {
        _cancellation?.Cancel();
        if (_runner is not null)
        {
            await _docker.StopContainerAsync(_runner.SpritzContainerName);
        }
    }

    /// <summary>Appends a line to the information pane and to the workflow log on disk.</summary>
    private void AppendInformation(string line)
    {
        string scrubbed = AnsiEscape.Replace(line ?? string.Empty, string.Empty);
        Information += scrubbed + Environment.NewLine;

        if (_runner is not null && !string.IsNullOrEmpty(_runner.PathToWorkflow))
        {
            try
            {
                File.AppendAllText(_runner.PathToWorkflow, scrubbed + Environment.NewLine);
            }
            catch (IOException)
            {
                // The pane is the primary output; a locked log file must not stop the run.
            }
        }
    }
}
