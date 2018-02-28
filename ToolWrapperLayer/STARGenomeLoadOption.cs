namespace ToolWrapperLayer
{
    /// <summary>
    /// Options for loading STAR genome indices.
    /// </summary>
    public enum STARGenomeLoadOption
    {
        NoSharedMemory,
        LoadAndExit,
        LoadAndKeep,
        LoadAndRemove,
        Remove,
    }
}