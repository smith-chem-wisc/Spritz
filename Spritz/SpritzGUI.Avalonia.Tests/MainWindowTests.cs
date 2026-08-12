using System.Linq;
using Avalonia;
using Avalonia.Headless;
using Avalonia.Headless.NUnit;
using NUnit.Framework;
using Spritz.Avalonia;
using Spritz.Avalonia.Services;
using Spritz.Avalonia.ViewModels;
using Spritz.Avalonia.Views;

[assembly: AvaloniaTestApplication(typeof(SpritzTest.Avalonia.TestAppBuilder))]

namespace SpritzTest.Avalonia;

/// <summary>Boots the real App on Avalonia's headless platform, so no display is needed.</summary>
public static class TestAppBuilder
{
    public static AppBuilder BuildAvaloniaApp() => AppBuilder.Configure<App>()
        .UseHeadless(new AvaloniaHeadlessPlatformOptions());
}

/// <summary>
/// The WPF GUI cannot be tested at all - its logic lives in MainWindow.xaml.cs behind a Window. These
/// tests exist to show what the port buys: the window can be constructed without a display, and the
/// command construction that used to be inline is now reachable directly.
/// </summary>
public class MainWindowTests
{
    /// <summary>
    /// Constructs the real window with the real view model. This is what a plain build cannot prove:
    /// that the XAML resolves at runtime, including the DataGrid theme brought in by App.axaml.
    /// </summary>
    [AvaloniaTest]
    public void MainWindowBuildsWithoutADisplay()
    {
        var window = new MainWindow { DataContext = new MainWindowViewModel() };
        window.Show();

        Assert.That(window.IsVisible, Is.True);
        Assert.That(window.Title, Is.EqualTo("Spritz"));
    }
}

/// <summary>
/// Docker argument construction, which the WPF GUI builds as one shell string inside a click handler.
/// </summary>
public class DockerArgumentsTests
{
    [Test]
    public void PublishedImageIsPulledAndTaggedWithTheCompiledVersion()
    {
        Assert.That(DockerArguments.ShouldPull("smithlab/spritz"), Is.True);
        Assert.That(DockerArguments.ImageWithVersion("smithlab/spritz"),
            Is.EqualTo($"smithlab/spritz:{SpritzBackend.RunnerEngine.CurrentVersion}"));
    }

    [Test]
    public void ExplicitTagIsUsedVerbatim()
    {
        Assert.That(DockerArguments.ImageWithVersion("spritz:dev"), Is.EqualTo("spritz:dev"),
            "a name that already carries a tag must not have the version appended");
    }

    [Test]
    public void LocalImageIsNotPulled()
    {
        Assert.That(DockerArguments.ShouldPull("spritz:dev"), Is.False);
        Assert.That(DockerArguments.ShouldPull("spritz"), Is.False);
    }

    /// <summary>
    /// The heuristic the WPF code uses is dockerImageName.Contains("smithlab"), which also matches a
    /// local image called "smithlab-test" and would pull over it. This asserts the stricter prefix.
    /// </summary>
    [Test]
    public void NameThatMerelyContainsSmithlabIsNotTreatedAsPublished()
    {
        Assert.That(DockerArguments.ShouldPull("smithlab-test"), Is.False);
        Assert.That(DockerArguments.ShouldPull("my-smithlab/spritz"), Is.False);
    }

    /// <summary>
    /// The point of the port's process handling: paths are separate arguments, so a space in one
    /// cannot split it. The WPF version interpolates paths into a shell string inside triple quotes.
    /// </summary>
    [Test]
    public void PathsWithSpacesSurviveAsSingleArguments()
    {
        var arguments = DockerArguments.Run(
            "smithlab/spritz", "spritz-test",
            "/home/a b/analysis", "/home/a b/resources",
            new[] { "dotnet", "SpritzCMD.dll" });

        Assert.That(arguments, Has.Member("/home/a b/analysis:/app/spritz/results/"));
        Assert.That(arguments, Has.Member("/home/a b/resources:/app/spritz/resources"));
        Assert.That(arguments.Count(a => a == "-v"), Is.EqualTo(2));
        Assert.That(arguments.Any(a => a.Contains('"')), Is.False,
            "no argument should need quoting, because no shell ever parses these");
    }

    [Test]
    public void RunNamesTheContainerSoItCanBeStopped()
    {
        var arguments = DockerArguments.Run(
            "smithlab/spritz", "spritz-abc", "/analysis", "/resources", new[] { "dotnet" });

        int nameIndex = arguments.ToList().IndexOf("--name");
        Assert.That(nameIndex, Is.GreaterThanOrEqualTo(0));
        Assert.That(arguments[nameIndex + 1], Is.EqualTo("spritz-abc"));
    }
}

public class DockerSystemInfoTests
{
    [Test]
    public void TotalMemoryIsReadFromDockerSystemInfo()
    {
        const string info = "Server:\n Containers: 0\n CPUs: 10\n Total Memory: 7.653GiB\n Name: docker-desktop\n";
        Assert.That(MainWindowViewModel.ParseTotalMemoryGb(info), Is.EqualTo(7.653 * 1.07374).Within(0.001));
    }

    [Test]
    public void MissingMemoryLineReportsZeroRatherThanThrowing()
    {
        Assert.That(MainWindowViewModel.ParseTotalMemoryGb("Cannot connect to the Docker daemon"), Is.EqualTo(0));
    }
}
