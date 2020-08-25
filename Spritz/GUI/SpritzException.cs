using System;

namespace Spritz
{
    [Serializable]
    public class SpritzException : Exception
    {
        public SpritzException(string message) : base(message)
        {
        }
    }
}