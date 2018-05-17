using System.Collections.Generic;
using System.IO;

namespace ToolWrapperLayer
{
    /// <summary>
    /// BLAST+ is a set of local alignment tools for nucleotide (blastn, etc.) and protein (blastp, etc.) sequences.
    ///
    /// BLASTN USAGE (Nucleotide-Nucleotide BLAST 2.2.28+)
    ///  blastn[-h][-help][-import_search_strategy filename]
    ///    [-export_search_strategy filename][-task task_name][-db database_name]
    ///    [-dbsize num_letters][-gilist filename][-seqidlist filename]
    ///    [-negative_gilist filename][-entrez_query entrez_query]
    ///    [-db_soft_mask filtering_algorithm][-db_hard_mask filtering_algorithm]
    ///    [-subject subject_input_file][-subject_loc range][-query input_file]
    ///    [-out output_file][-evalue evalue][-word_size int_value]
    ///    [-gapopen open_penalty][-gapextend extend_penalty]
    ///    [-perc_identity float_value][-xdrop_ungap float_value]
    ///    [-xdrop_gap float_value][-xdrop_gap_final float_value]
    ///    [-searchsp int_value][-max_hsps_per_subject int_value][-penalty penalty]
    ///    [-reward reward][-no_greedy][-min_raw_gapped_score int_value]
    ///    [-template_type type][-template_length int_value][-dust DUST_options]
    ///    [-filtering_db filtering_database]
    ///    [-window_masker_taxid window_masker_taxid]
    ///    [-window_masker_db window_masker_db][-soft_masking soft_masking]
    ///    [-ungapped][-culling_limit int_value][-best_hit_overhang float_value]
    ///    [-best_hit_score_edge float_value][-window_size int_value]
    ///    [-off_diagonal_range int_value][-use_index boolean][-index_name string]
    ///    [-lcase_masking][-query_loc range][-strand strand][-parse_deflines]
    ///    [-outfmt format][-show_gis][-num_descriptions int_value]
    ///    [-num_alignments int_value][-html][-max_target_seqs num_sequences]
    ///    [-num_threads int_value][-remote][-version]
    ///
    /// BLASTP USAGE (Protein-Protein BLAST 2.2.28+)
    ///  blastp[-h][-help][-import_search_strategy filename]
    ///    [-export_search_strategy filename][-task task_name][-db database_name]
    ///    [-dbsize num_letters][-gilist filename][-seqidlist filename]
    ///    [-negative_gilist filename][-entrez_query entrez_query]
    ///    [-db_soft_mask filtering_algorithm][-db_hard_mask filtering_algorithm]
    ///    [-subject subject_input_file][-subject_loc range][-query input_file]
    ///    [-out output_file][-evalue evalue][-word_size int_value]
    ///    [-gapopen open_penalty][-gapextend extend_penalty]
    ///    [-xdrop_ungap float_value][-xdrop_gap float_value]
    ///    [-xdrop_gap_final float_value][-searchsp int_value]
    ///    [-max_hsps_per_subject int_value][-seg SEG_options]
    ///    [-soft_masking soft_masking][-matrix matrix_name]
    ///    [-threshold float_value][-culling_limit int_value]
    ///    [-best_hit_overhang float_value][-best_hit_score_edge float_value]
    ///    [-window_size int_value][-lcase_masking][-query_loc range]
    ///    [-parse_deflines][-outfmt format][-show_gis]
    ///    [-num_descriptions int_value][-num_alignments int_value][-html]
    ///    [-max_target_seqs num_sequences][-num_threads int_value][-ungapped]
    ///    [-remote][-comp_based_stats compo][-use_sw_tback][-version]
    /// </summary>
    public class BLASTWrapper :
        IInstallable
    {
        #region Installation Methods

        /// <summary>
        /// Writes an install script for BLAST+.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteInstallScript(string spritzDirectory)
        {
            string scriptPath = WrapperUtility.GetInstallationScriptPath(spritzDirectory, "InstallBLAST.bash");
            WrapperUtility.GenerateScript(scriptPath, new List<string>
            {
                WrapperUtility.ChangeToToolsDirectoryCommand(spritzDirectory),
                "if [ ! -d ncbi-blast-2.7.1+ ]; then",
                "  wget --no-check ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz",
                "  tar -xvf ncbi-blast-2.7.1+-x64-linux.tar.gz; rm ncbi-blast-2.7.1+-x64-linux.tar.gz",
                "  cp ncbi-blast-2.7.1+/bin/* /usr/local/bin",
                "fi"
            });
            return scriptPath;
        }

        /// <summary>
        /// Writes a script for removing bedtools.
        /// </summary>
        /// <param name="spritzDirectory"></param>
        /// <returns></returns>
        public string WriteRemoveScript(string spritzDirectory)
        {
            return null;
        }

        #endregion Installation Methods
    }
}