using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using WorkflowLayer;

namespace SpritzGUI
{
    /// <summary>
    /// Interaction logic for STARAlignWorkFlowWindows.xaml
    /// </summary>
    public partial class STARAlignWorkFlowWindows : Window
    {
        public STARAlignWorkFlowWindows()
        {
            InitializeComponent();
        }

        internal STARAlignmentFlow TheTask { get; private set; }

        protected void cancelButton_Click(object sender, RoutedEventArgs e)
        {

        }

        protected void saveButton_Click(object sender, RoutedEventArgs e)
        {

        }
    }
}
