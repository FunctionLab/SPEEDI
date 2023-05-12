#' Perform differential expression analysis (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param metadata_df Dataframe containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A list containing differential expression analyses
#' @examples
#' \dontrun{differential_expression_results <- RunDE_RNA(sc_obj = sc_obj, metadata_df = metadata_df)}
#' @export
RunDE_RNA <- function(sc_obj, metadata_df, output_dir = getwd(), log_flag = FALSE) {
  species <- tolower(species)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Running differential expression analysis (RNA)", log_flag)
  Seurat::DefaultAssay(sc_obj) <- "RNA"
  for(metadata_attribute in colnames(metadata_df)) {
    current_de <- Libra::run_de(sc_obj, replicate_col = "sample",
                         cell_type_col = "predicted_celltype_majority_vote", label_col = metadata_attribute)
    utils::write.table(current_de, file = paste0(output_dir, metadata_attribute, ".tsv"), sep = "\t", quote = FALSE)
  }
  print_SPEEDI("Differential expression analysis complete", log_flag)
  gc()
  return(TRUE)
}
