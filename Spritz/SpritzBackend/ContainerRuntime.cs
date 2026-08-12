namespace SpritzBackend
{
    /// <summary>
    /// Which container runtime Spritz drives. All three run the same OCI image, so this only changes
    /// how the command is spelled.
    /// </summary>
    public enum ContainerRuntime
    {
        /// <summary>Default. Apache-2.0, daemonless, and installable without Docker Desktop.</summary>
        Podman,

        /// <summary>Kept because it is what existing users already have installed.</summary>
        Docker,

        /// <summary>For HPC, where a root daemon is usually not permitted and Apptainer is.</summary>
        Apptainer,
    }

    public static class ContainerRuntimes
    {
        /// <summary>Parses a runtime name, falling back to the default rather than throwing.</summary>
        public static ContainerRuntime Parse(string name)
        {
            if (string.IsNullOrWhiteSpace(name))
                return ContainerRuntime.Podman;

            switch (name.Trim().ToLowerInvariant())
            {
                case "docker": return ContainerRuntime.Docker;
                case "apptainer":
                case "singularity": return ContainerRuntime.Apptainer;
                default: return ContainerRuntime.Podman;
            }
        }

        /// <summary>The executable name, which is also the CLI Spritz has to speak.</summary>
        public static string Executable(this ContainerRuntime runtime) =>
            runtime == ContainerRuntime.Apptainer ? "apptainer"
            : runtime == ContainerRuntime.Docker ? "docker"
            : "podman";

        /// <summary>
        /// Docker and Podman take the same flags for everything Spritz uses, so they share a command
        /// shape. Apptainer does not: it has no daemon, so nothing to name or stop, it binds rather
        /// than mounts, and it reads OCI images through a docker:// URI.
        /// </summary>
        public static bool IsDockerCompatible(this ContainerRuntime runtime) =>
            runtime != ContainerRuntime.Apptainer;
    }
}
