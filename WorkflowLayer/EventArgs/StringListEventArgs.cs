using System;
using System.Collections.Generic;

namespace WorkflowLayer
{
    public class StringListEventArgs : EventArgs
    {
        #region Public Constructors

        public StringListEventArgs(List<string> stringList)
        {
            this.StringList = stringList;
        }

        #endregion Public Constructors

        #region Public Properties

        public IEnumerable<string> StringList { get; }

        #endregion Public Properties
    }
}