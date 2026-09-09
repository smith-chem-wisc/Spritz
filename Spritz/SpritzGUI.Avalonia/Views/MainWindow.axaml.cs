using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Avalonia.Controls;
using Avalonia.Interactivity;
using Avalonia.Markup.Xaml;
using Avalonia.Platform.Storage;
using Spritz.Avalonia.ViewModels;

namespace Spritz.Avalonia.Views;

/// <summary>
/// The view holds only what genuinely needs a window: file pickers, opening a browser, and the
/// workflow dialog. Everything else lives in the view model.
/// </summary>
public partial class MainWindow : Window
{
    private MainWindowViewModel ViewModel => (MainWindowViewModel)DataContext;

    public MainWindow() => AvaloniaXamlLoader.Load(this);

    /// <summary>
    /// Avalonia's StorageProvider replaces Microsoft.Win32.OpenFileDialog, which is Windows-only.
    /// </summary>
    private async void OnAddFastqClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFile> files = await StorageProvider.OpenFilePickerAsync(new FilePickerOpenOptions
        {
            Title = "Select FASTQ files",
            AllowMultiple = true,
            FileTypeFilter = new[]
            {
                new FilePickerFileType("FASTQ") { Patterns = new[] { "*.fastq", "*.fastq.gz" } },
            },
        });

        ViewModel.AddFastqFiles(files.Select(f => f.Path.LocalPath));
    }

    private async void OnBrowseOutputClick(object sender, RoutedEventArgs e)
    {
        IReadOnlyList<IStorageFolder> folders = await StorageProvider.OpenFolderPickerAsync(
            new FolderPickerOpenOptions { Title = "Select the analysis directory", AllowMultiple = false });

        if (folders.Count > 0)
        {
            ViewModel.OutputFolder = folders[0].Path.LocalPath;
        }
    }

    private async void OnCreateWorkflowClick(object sender, RoutedEventArgs e)
    {
        // The workflow dialog is the remaining piece of the port. Until it exists, say so plainly
        // rather than appearing to work.
        await ShowMessageAsync("Not yet ported",
            "The workflow options dialog has not been ported yet. Until it is, build a workflow in the "
            + "Windows GUI, or run the container directly with SpritzCMD.");
    }

    private void OnWikiClick(object sender, RoutedEventArgs e) =>
        OpenUrl("https://github.com/smith-chem-wisc/Spritz/wiki");

    private void OnContactClick(object sender, RoutedEventArgs e) =>
        OpenUrl("https://github.com/smith-chem-wisc/Spritz/issues");

    /// <summary>Process.Start with UseShellExecute opens the platform's default browser on all three OSes.</summary>
    private static void OpenUrl(string url)
    {
        try
        {
            Process.Start(new ProcessStartInfo(url) { UseShellExecute = true });
        }
        catch (Exception)
        {
            // Opening a browser is a convenience; failing to do so must not take the window down.
        }
    }

    /// <summary>A minimal replacement for WPF's MessageBox, which Avalonia does not provide.</summary>
    private async System.Threading.Tasks.Task ShowMessageAsync(string title, string message)
    {
        var dialog = new Window
        {
            Title = title,
            Width = 420,
            SizeToContent = SizeToContent.Height,
            WindowStartupLocation = WindowStartupLocation.CenterOwner,
        };
        var ok = new Button { Content = "OK", HorizontalAlignment = global::Avalonia.Layout.HorizontalAlignment.Right };
        ok.Click += (_, _) => dialog.Close();
        dialog.Content = new StackPanel
        {
            Margin = new global::Avalonia.Thickness(16),
            Spacing = 12,
            Children =
            {
                new TextBlock { Text = message, TextWrapping = global::Avalonia.Media.TextWrapping.Wrap },
                ok,
            },
        };
        await dialog.ShowDialog(this);
    }
}
