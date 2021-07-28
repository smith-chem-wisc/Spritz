using System;

namespace SpritzBackend
{
    [Serializable]
    public class SpritzException : Exception
    {
        public SpritzException(string message) : base(message)
        {
        }
    }
}