namespace GeneProtGUI
{
    partial class ProteoformDatabaseEngine
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.btn_writeXmlFromXml = new System.Windows.Forms.Button();
            this.btn_writeXmlFromFasta = new System.Windows.Forms.Button();
            this.tb_input = new System.Windows.Forms.TextBox();
            this.tb_output = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // btn_writeXmlFromXml
            // 
            this.btn_writeXmlFromXml.Location = new System.Drawing.Point(221, 256);
            this.btn_writeXmlFromXml.Name = "btn_writeXmlFromXml";
            this.btn_writeXmlFromXml.Size = new System.Drawing.Size(216, 23);
            this.btn_writeXmlFromXml.TabIndex = 0;
            this.btn_writeXmlFromXml.Text = "Write xml from xml";
            this.btn_writeXmlFromXml.UseVisualStyleBackColor = true;
            // 
            // btn_writeXmlFromFasta
            // 
            this.btn_writeXmlFromFasta.Location = new System.Drawing.Point(221, 285);
            this.btn_writeXmlFromFasta.Name = "btn_writeXmlFromFasta";
            this.btn_writeXmlFromFasta.Size = new System.Drawing.Size(216, 23);
            this.btn_writeXmlFromFasta.TabIndex = 1;
            this.btn_writeXmlFromFasta.Text = "Write xml from fasta";
            this.btn_writeXmlFromFasta.UseVisualStyleBackColor = true;
            // 
            // tb_input
            // 
            this.tb_input.Location = new System.Drawing.Point(78, 104);
            this.tb_input.Name = "tb_input";
            this.tb_input.Size = new System.Drawing.Size(100, 20);
            this.tb_input.TabIndex = 2;
            // 
            // tb_output
            // 
            this.tb_output.Location = new System.Drawing.Point(78, 139);
            this.tb_output.Name = "tb_output";
            this.tb_output.Size = new System.Drawing.Size(100, 20);
            this.tb_output.TabIndex = 3;
            // 
            // ProteoformDatabaseEngine
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(629, 463);
            this.Controls.Add(this.tb_output);
            this.Controls.Add(this.tb_input);
            this.Controls.Add(this.btn_writeXmlFromFasta);
            this.Controls.Add(this.btn_writeXmlFromXml);
            this.Name = "ProteoformDatabaseEngine";
            this.Text = "ProteoformDatabaseEngine";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button btn_writeXmlFromXml;
        private System.Windows.Forms.Button btn_writeXmlFromFasta;
        private System.Windows.Forms.TextBox tb_input;
        private System.Windows.Forms.TextBox tb_output;
    }
}

