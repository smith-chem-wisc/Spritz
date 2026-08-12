using Avalonia;

namespace Spritz.Avalonia;

internal static class Program
{
    // Avalonia needs a plain Main; nothing here is platform specific.
    [System.STAThread]
    public static void Main(string[] args) => BuildAvaloniaApp()
        .StartWithClassicDesktopLifetime(args);

    public static AppBuilder BuildAvaloniaApp() => AppBuilder.Configure<App>()
        .UsePlatformDetect()
        .LogToTrace();
}
