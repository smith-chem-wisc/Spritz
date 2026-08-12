using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Spritz.Avalonia.ViewModels;

namespace SpritzTest.Avalonia;

/// <summary>
/// Behaviour tests modelled on the patterns in MetaMorpheus's Test/GuiTests suite, which the first
/// pass here missed. Those 299 tests consistently cover four things this file now covers too:
/// that property changes raise notifications, that commands guard their preconditions, that
/// constructors establish sensible defaults, and that bad input is rejected rather than absorbed.
///
/// Without the notification tests in particular, a binding can silently stop updating and nothing
/// fails - the window just stops reflecting reality.
/// </summary>
public class ViewModelNotificationTests
{
    /// <summary>Collects the property names a view model reports as changed.</summary>
    private static List<string> RecordChanges(INotifyPropertyChanged source)
    {
        var changed = new List<string>();
        source.PropertyChanged += (_, e) => changed.Add(e.PropertyName);
        return changed;
    }

    [Test]
    public void EditableStateRaisesPropertyChanged()
    {
        var viewModel = new MainWindowViewModel();
        var changed = RecordChanges(viewModel);

        viewModel.SraAccession = "SRR13737862";
        viewModel.OutputFolder = "/tmp/analysis";
        viewModel.DockerImage = "spritz:dev";

        Assert.That(changed, Does.Contain(nameof(MainWindowViewModel.SraAccession)));
        Assert.That(changed, Does.Contain(nameof(MainWindowViewModel.OutputFolder)));
        Assert.That(changed, Does.Contain(nameof(MainWindowViewModel.DockerImage)));
    }

    [Test]
    public void SettingAPropertyToItsCurrentValueRaisesNothing()
    {
        var viewModel = new MainWindowViewModel { OutputFolder = "/tmp/analysis" };
        var changed = RecordChanges(viewModel);

        viewModel.OutputFolder = "/tmp/analysis";

        Assert.That(changed, Is.Empty, "an unchanged value should not notify, or the UI churns");
    }

    [Test]
    public void RunStateRaisesPropertyChangedSoButtonsCanEnableAndDisable()
    {
        var viewModel = new MainWindowViewModel();
        var changed = RecordChanges(viewModel);

        viewModel.IsRunning = true;

        Assert.That(changed, Does.Contain(nameof(MainWindowViewModel.IsRunning)),
            "the Run and Cancel buttons bind to IsRunning, so it has to notify");
    }
}

public class SraEntryTests
{
    [TestCase("SRR13737862")]
    [TestCase("ERR1234567")]
    [TestCase("srr13737862")]
    public void RecognisableAccessionsAreAccepted(string accession)
    {
        var viewModel = new MainWindowViewModel { SraAccession = accession };
        viewModel.AddPairedEndSraCommand.Execute(null);

        Assert.That(viewModel.Sras.Select(s => s.Name), Does.Contain(accession));
    }

    [TestCase("not-an-accession")]
    [TestCase("12345")]
    [TestCase("")]
    public void UnrecognisableAccessionsAreRejectedWithAnExplanation(string accession)
    {
        var viewModel = new MainWindowViewModel { SraAccession = accession };
        viewModel.AddPairedEndSraCommand.Execute(null);

        Assert.That(viewModel.Sras, Is.Empty);
        if (accession.Length > 0)
        {
            Assert.That(viewModel.Information, Does.Contain(accession),
                "a rejected accession should be named, so the user knows which one was wrong");
        }
    }

    [Test]
    public void PairedAndSingleEndAreRecordedDistinctly()
    {
        var viewModel = new MainWindowViewModel { SraAccession = "SRR111" };
        viewModel.AddPairedEndSraCommand.Execute(null);
        viewModel.SraAccession = "SRR222";
        viewModel.AddSingleEndSraCommand.Execute(null);

        Assert.That(viewModel.Sras.Single(s => s.Name == "SRR111").IsPairedEnd, Is.True);
        Assert.That(viewModel.Sras.Single(s => s.Name == "SRR222").IsPairedEnd, Is.False);
    }

    [Test]
    public void TheSameAccessionIsNotAddedTwice()
    {
        var viewModel = new MainWindowViewModel { SraAccession = "SRR13737862" };
        viewModel.AddPairedEndSraCommand.Execute(null);
        viewModel.SraAccession = "SRR13737862";
        viewModel.AddPairedEndSraCommand.Execute(null);

        Assert.That(viewModel.Sras, Has.Count.EqualTo(1),
            "adding a duplicate accession would make the pipeline download it twice");
    }

    [Test]
    public void AddingClearsTheEntryBoxSoTheNextOneCanBeTyped()
    {
        var viewModel = new MainWindowViewModel { SraAccession = "SRR13737862" };
        viewModel.AddPairedEndSraCommand.Execute(null);

        Assert.That(viewModel.SraAccession, Is.Empty);
    }

    [Test]
    public void ClearRemovesEverything()
    {
        var viewModel = new MainWindowViewModel { SraAccession = "SRR1" };
        viewModel.AddPairedEndSraCommand.Execute(null);

        viewModel.ClearSrasCommand.Execute(null);

        Assert.That(viewModel.Sras, Is.Empty);
    }
}

public class FastqEntryTests
{
    /// <summary>
    /// Mate-pair detection reads the _1 and _2 suffixes off the filename. Getting it wrong silently
    /// mispairs the reads, which is the sort of thing that produces a plausible but wrong result, so
    /// it is worth pinning even though it looks trivial.
    /// </summary>
    [Test]
    public void MatePairIsReadFromTheFilenameSuffix()
    {
        var viewModel = new MainWindowViewModel();
        viewModel.AddFastqFiles(new[] { "/data/TK12_1.fastq", "/data/TK12_2.fastq" });

        Assert.That(viewModel.Fastqs.Single(f => f.FileName == "TK12_1").MatePair, Is.EqualTo("1"));
        Assert.That(viewModel.Fastqs.Single(f => f.FileName == "TK12_2").MatePair, Is.EqualTo("2"));
    }

    [Test]
    public void GzippedFastqsKeepTheirLogicalName()
    {
        var viewModel = new MainWindowViewModel();
        viewModel.AddFastqFiles(new[] { "/data/sample_1.fastq.gz" });

        Assert.That(viewModel.Fastqs, Has.Count.EqualTo(1));
        Assert.That(viewModel.Fastqs[0].FileName, Is.EqualTo("sample_1"),
            "both extensions should come off, otherwise the name carries .fastq");
    }

    [TestCase("/data/notes.txt")]
    [TestCase("/data/reads.bam")]
    public void NonFastqFilesAreIgnored(string path)
    {
        var viewModel = new MainWindowViewModel();
        viewModel.AddFastqFiles(new[] { path });

        Assert.That(viewModel.Fastqs, Is.Empty);
    }

    [Test]
    public void TheSameFileIsNotAddedTwice()
    {
        var viewModel = new MainWindowViewModel();
        viewModel.AddFastqFiles(new[] { "/data/sample_1.fastq" });
        viewModel.AddFastqFiles(new[] { "/data/sample_1.fastq" });

        Assert.That(viewModel.Fastqs, Has.Count.EqualTo(1));
    }

    /// <summary>
    /// The analysis directory defaults to where the FASTQs are, because the pipeline expects them
    /// under it - #237's reporter hit trouble partly through a mismatch here.
    /// </summary>
    [Test]
    public void TheOutputFolderDefaultsToWhereTheFastqsAre()
    {
        var viewModel = new MainWindowViewModel();
        viewModel.AddFastqFiles(new[] { Path.Combine("/data", "reads", "sample_1.fastq") });

        Assert.That(viewModel.OutputFolder, Is.EqualTo(Path.Combine("/data", "reads")));
    }

    [Test]
    public void AnExistingOutputFolderIsNotOverwritten()
    {
        var viewModel = new MainWindowViewModel { OutputFolder = "/chosen/by/the/user" };
        viewModel.AddFastqFiles(new[] { "/data/sample_1.fastq" });

        Assert.That(viewModel.OutputFolder, Is.EqualTo("/chosen/by/the/user"),
            "adding files must not silently move an analysis directory the user picked");
    }
}

public class RunPreconditionTests
{
    [Test]
    public void RunIsRefusedWithNoWorkflow()
    {
        var viewModel = new MainWindowViewModel { OutputFolder = "/tmp/out" };

        Assert.That(viewModel.CanRun, Is.False, "there is no workflow to run yet");
        viewModel.RunWorkflowCommand.Execute(null);
        Assert.That(viewModel.Information, Does.Contain("workflow"),
            "refusing should say why, not fail silently");
    }

    [Test]
    public void RunIsRefusedWithNoAnalysisDirectory()
    {
        var viewModel = new MainWindowViewModel();

        Assert.That(viewModel.CanRun, Is.False);
    }

    [Test]
    public void ADefaultDockerImageIsSetSoTheGuiWorksUnconfigured()
    {
        Assert.That(new MainWindowViewModel().DockerImage, Is.EqualTo("smithlab/spritz"));
    }

    [Test]
    public void DockerStatusStartsAsAnInProgressMessageRatherThanBlank()
    {
        Assert.That(new MainWindowViewModel().DockerStatus, Is.Not.Empty,
            "a blank status bar at startup looks like a broken window");
    }
}
