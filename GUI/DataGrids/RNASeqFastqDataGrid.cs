using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SpritzGUI
{
    class RNASeqFastqDataGrid
    {
        public RNASeqFastqDataGrid(string filePath)
        {
            Use = true;
            FilePath = filePath;
        }

        public bool Use { get; set; }
        public bool InProgress { get; private set; }
        public string FilePath { get; set; }

        public void SetInProgress(bool inProgress)
        {
            InProgress = inProgress;
        }
    }
}
