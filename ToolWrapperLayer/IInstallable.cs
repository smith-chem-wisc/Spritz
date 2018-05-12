namespace ToolWrapperLayer
{
    public interface IInstallable
    {
        /// <summary>
        /// Writes an installation script
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns>Path of installation script (Windows-formatted)</returns>
        string WriteInstallScript(string spritzDirectory);

        /// <summary>
        /// Writes a cleaning script
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns>Path of cleaning script (Windows-formatted)</returns>
        string WriteRemoveScript(string spritzDirectory);
    }
}
