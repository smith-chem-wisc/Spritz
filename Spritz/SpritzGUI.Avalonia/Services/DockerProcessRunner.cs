using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using System.Threading.Tasks;

namespace Spritz.Avalonia.Services;

/// <summary>
/// Runs docker directly, rather than handing a command string to a shell.
///
/// The WPF GUI sets proc.StartInfo.FileName = "Powershell.exe" and passes the whole
/// `docker pull ...; docker run ...; docker stop ...` string as arguments. That only works on
/// Windows, and it relies on a shell to split and unquote the arguments - which is also why the
/// existing command interpolates volume paths inside triple quotes. Passing an argument list
/// instead removes both the platform dependency and the quoting question: no shell ever sees the
/// path, so spaces and quotes in it are not special.
/// </summary>
public sealed class DockerProcessRunner
{
    private readonly Action<string> _onOutput;

    public DockerProcessRunner(Action<string> onOutput) => _onOutput = onOutput;

    /// <summary>Whether a docker binary is present and responding, and what it said.</summary>
    public async Task<(bool Available, string Info)> QuerySystemInfoAsync(CancellationToken token = default)
    {
        try
        {
            var (exitCode, output) = await CaptureAsync("docker", new[] { "system", "info" }, token);
            bool available = exitCode == 0 && !string.IsNullOrWhiteSpace(output);
            return (available, output);
        }
        catch (Exception exception)
        {
            // A missing binary throws rather than returning a code, and that is the common case for
            // a user who has not installed Docker, so it is reported rather than thrown onward.
            return (false, exception.Message);
        }
    }

    /// <summary>Streams `docker run` output line by line until the container exits.</summary>
    public Task<int> RunAsync(IReadOnlyList<string> arguments, CancellationToken token = default) =>
        StreamAsync("docker", arguments, token);

    /// <summary>Stops a named container. Used on cancel and on window close.</summary>
    public async Task StopContainerAsync(string containerName, CancellationToken token = default)
    {
        if (string.IsNullOrWhiteSpace(containerName))
        {
            return;
        }
        await CaptureAsync("docker", new[] { "stop", containerName }, token);
    }

    private async Task<int> StreamAsync(string fileName, IReadOnlyList<string> arguments, CancellationToken token)
    {
        using var process = new Process { StartInfo = StartInfo(fileName, arguments) };
        process.OutputDataReceived += (_, e) => Emit(e.Data);
        process.ErrorDataReceived += (_, e) => Emit(e.Data);

        process.Start();
        process.BeginOutputReadLine();
        process.BeginErrorReadLine();

        await using (token.Register(() => TryKill(process)))
        {
            await process.WaitForExitAsync(CancellationToken.None);
        }
        return process.ExitCode;
    }

    private static async Task<(int ExitCode, string Output)> CaptureAsync(
        string fileName, IReadOnlyList<string> arguments, CancellationToken token)
    {
        using var process = new Process { StartInfo = StartInfo(fileName, arguments) };
        process.Start();
        string output = await process.StandardOutput.ReadToEndAsync(token);
        string error = await process.StandardError.ReadToEndAsync(token);
        await process.WaitForExitAsync(token);
        return (process.ExitCode, string.IsNullOrWhiteSpace(output) ? error : output);
    }

    private static ProcessStartInfo StartInfo(string fileName, IReadOnlyList<string> arguments)
    {
        var info = new ProcessStartInfo
        {
            FileName = fileName,
            UseShellExecute = false,
            RedirectStandardOutput = true,
            RedirectStandardError = true,
            CreateNoWindow = true,
        };
        foreach (string argument in arguments)
        {
            info.ArgumentList.Add(argument);
        }
        return info;
    }

    private void Emit(string line)
    {
        if (!string.IsNullOrEmpty(line))
        {
            _onOutput(line);
        }
    }

    private static void TryKill(Process process)
    {
        try
        {
            if (!process.HasExited)
            {
                process.Kill(entireProcessTree: true);
            }
        }
        catch (InvalidOperationException)
        {
            // Already gone between the check and the kill; nothing to do.
        }
    }
}
