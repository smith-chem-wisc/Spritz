using System.Threading.Tasks;
using ToolWrapperLayer;

namespace WorkflowLayer
{
    public class InstallFlow
    {
        public static void Run(string binDirectory)
        {
            WrapperUtility.Install(binDirectory);
            Parallel.Invoke
            (
                () => GATKWrapper.Install(binDirectory),
                () => RSeQCWrapper.Install(binDirectory),
                () => ScalpelWrapper.Install(binDirectory),
                () => SkewerWrapper.Install(binDirectory),
                () => SlnckyWrapper.Install(binDirectory),
                () => SRAToolkitWrapper.Install(binDirectory),
                () => STARWrapper.Install(binDirectory),
                () => STARFusionWrapper.Install(binDirectory)
            );
        }
    }
}
