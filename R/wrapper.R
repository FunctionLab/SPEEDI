#' Wrapper function for running the SPEEDI pipeline
#'
#' @param tissue tissue of data (used to map cell types via reference)
#' @param data_path path to where data is located
#' @param sample_id_list list of sample names (optional - if not provided, will select all samples found recursively in data_path)
#' @param human flag to indicate whether we're processing human or mouse data
#' @param record_doublets flag to indicate whether we're recording doublets in the data (using scDblFinder)
#' @return A Seurat object that has been processed through the SPEEDI pipeline
#' @export
run_SPEEDI <- function(tissue, data_path = getwd(), sample_id_list = NULL, human = TRUE, record_doublets = FALSE) {
  all_sc_exp_matrices <- Read_h5(data_path, sample_id_list)
  sc_obj <- FilterRawData(all_sc_exp_matrices, human, record_doublets)
  rm(all_sc_exp_matrices)
  sc_obj <- InitialProcessing(sc_obj, human)
  sc_obj <- InferBatches(sc_obj)
  sc_obj <- IntegrateByBatch(sc_obj)
  sc_obj <- VisualizeIntegration(sc_obj)
  reference <- LoadReference(tissue, human)
  sc_obj <- MapCellTypes(sc_obj, reference)
  return(sc_obj)
}
