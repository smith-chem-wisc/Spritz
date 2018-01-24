namespace ToolWrapperLayer
{
    public interface IInstallable
    {
        string WriteInstallScript(string binDirectory);
        string WriteRemoveScript(string binDirectory);
    }
}
