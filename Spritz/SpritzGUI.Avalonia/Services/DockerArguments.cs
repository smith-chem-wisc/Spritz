using System;
using System.Collections.Generic;
using SpritzBackend;

namespace Spritz.Avalonia.Services;

/// <summary>
/// Builds the docker argument list, as a list rather than a shell string.
///
/// This is the argument-list counterpart of RunnerEngine.GenerateCommandsDry, which returns
/// `docker pull ...; docker run ...; docker stop ...` as one string for a shell to split. Keeping it
/// here rather than in SpritzBackend for now, so the draft does not change behaviour the Windows GUI
/// depends on; it belongs in SpritzBackend once both front ends can use it.
/// </summary>
public static class DockerArguments
{
    /// <summary>Whether this image name refers to the published image and so should be pulled first.</summary>
    public static bool ShouldPull(string dockerImageName) =>
        dockerImageName.StartsWith("smithlab/", StringComparison.Ordinal);

    /// <summary>The image reference, with the compiled version appended when no tag was given.</summary>
    public static string ImageWithVersion(string dockerImageName) =>
        dockerImageName.Contains(':') ? dockerImageName : $"{dockerImageName}:{RunnerEngine.CurrentVersion}";

    public static IReadOnlyList<string> Pull(string dockerImageName) =>
        new[] { "pull", ImageWithVersion(dockerImageName) };

    /// <summary>
    /// `docker run` with the analysis and resources directories mounted. The paths go in as
    /// individual arguments, so a space or a quote in them needs no escaping and cannot split.
    /// </summary>
    public static IReadOnlyList<string> Run(
        string dockerImageName, string containerName, string analysisDirectory,
        string resourcesDirectory, IEnumerable<string> containerCommand)
    {
        var arguments = new List<string>
        {
            "run", "--rm", "--user=root",
            "--name", containerName,
            "-v", $"{analysisDirectory}:/app/spritz/results/",
            "-v", $"{resourcesDirectory}:/app/spritz/resources",
            ImageWithVersion(dockerImageName),
        };
        arguments.AddRange(containerCommand);
        return arguments;
    }

    /// <summary>The in-container command: conda run, then SpritzCMD with the user's options.</summary>
    public static IReadOnlyList<string> SpritzCmd(SpritzOptions options)
    {
        var arguments = new List<string>
        {
            "conda", "run", "--no-capture-output", "--live-stream", "dotnet", "SpritzCMD.dll",
        };
        // SpritzOptionStrings builds a single argument string; split on spaces it does not own.
        foreach (string argument in SpritzOptionStrings.GenerateSpritzCMDArgs(options)
                     .Split(' ', StringSplitOptions.RemoveEmptyEntries))
        {
            arguments.Add(argument.Replace("\"\"\"", string.Empty));
        }
        return arguments;
    }

    /// <summary>A readable rendering of a command, for the information pane only.</summary>
    public static string Describe(IEnumerable<string> arguments) => "docker " + string.Join(" ", arguments);
}
