#' Wrapper function for running the SPEEDI pipeline
#'
#' @param tissue tissue of data (used to map cell types via reference)
#' @param data_path path to where data is located
#' @param reference_dir path to dir where reference is located
#' @param reference_file_name name of reference file
#' @param output_dir path to dir where output will be saved
#' @param sample_id_list list of sample names (optional - if not provided, will select all samples found recursively in data_path)
#' @param species flag to indicate whether we're processing human or mouse data
#' @param record_doublets flag to indicate whether we're recording doublets in the data (using scDblFinder)
#' @return A Seurat object that has been processed through the SPEEDI pipeline
#' @export
run_SPEEDI <- function(tissue, data_path = getwd(), reference_dir = getwd(), reference_file_name = NULL, output_dir = getwd(), sample_id_list = NULL, species = "human", record_doublets = FALSE) {
  # Create output_dir if it doesn't already exist
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  # Add "/" to end of output_dir if not already present
  last_char_of_output_dir_path <- substr(output_dir, nchar(output_dir), nchar(output_dir))
  if(last_char_of_output_dir_path != "/") {
    output_dir <- paste0(output_dir, "/")
  }
  # Create log file
  log_file_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
  log_file_name <- gsub(":", "-", log_file_name)
  log_file <- logr::log_open(paste0(output_dir, log_file_name), logdir = FALSE)
  # Stage 1 - read in data
  all_sc_exp_matrices <- Read_h5(data_path, sample_id_list, log_flag = TRUE)
  # Stage 2 - filter data
  sc_obj <- FilterRawData(all_sc_exp_matrices, species, record_doublets, log_flag = TRUE)
  rm(all_sc_exp_matrices)
  # Stage 3 - process data
  sc_obj <- InitialProcessing(sc_obj, species, log_flag = TRUE)
  # Stage 4 - find batches
  sc_obj <- InferBatches(sc_obj, log_flag = TRUE)
  # Stage 5 - integrate batches
  sc_obj <- IntegrateByBatch(sc_obj, log_flag = TRUE)
  # Stage 6 - create UMAP of integration (and prep for FindMarkers)
  sc_obj <- VisualizeIntegration(sc_obj, log_flag = TRUE)
  # Stage 7 - Load reference
  reference <- LoadReferenceSPEEDI(tissue, species, reference_dir, reference_file_name, log_flag = TRUE)
  # Stage 8 - Map cell types from reference onto query data
  sc_obj <- MapCellTypes(sc_obj, reference, log_flag = TRUE)
  return(sc_obj)
}
