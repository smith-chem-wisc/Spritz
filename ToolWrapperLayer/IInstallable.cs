namespace ToolWrapperLayer
{
    public interface IInstallable
    {
        string WriteInstallScript(string spritzDirectory);
        string WriteRemoveScript(string spritzDirectory);
    }
}
