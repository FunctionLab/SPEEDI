#' Find overlap between true multiome cells (RNA / ATAC)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param proj ArchR project associated with data
#' @param data_modality data modality to return (in Seurat object or ArchR project)
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object or ArchR project with only those cells that pass QC in the other assay
#' @examples
#' \dontrun{sc_obj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "RNA")}
#' \dontrun{proj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "ATAC")}
#' @export
FindMultiomeOverlap <- function(sc_obj, proj, data_modality = "RNA", log_flag = FALSE) {
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 9: Finding overlap between RNA and ATAC for multiome data", log_flag)
  print_SPEEDI(paste0("data_modality is: ", data_modality), log_flag)
  if(data_modality == "ATAC") {
    proj <- ArchR::subsetArchRProject(ArchRProj = proj,
                               cells = colnames(sc_obj),
                               force = TRUE)
  } else if(data_modality == "RNA") {
    sc_obj <- sc_obj[,ArchR::getCellNames(proj)]
  } else {
    print_SPEEDI("Error: you have selected an invalid data modality (you must choose RNA or ATAC)", log_flag)
    stop()
  }
  print_SPEEDI("Step 9: Complete", log_flag)
  gc()
  if(data_modality == "ATAC") {
    return(proj)
  } else {
    return(sc_obj)
  }
}

#' Transfer cell type labels from RNA to ATAC (only works for true multiome since same cells in both assays)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param proj ArchR project associated with data
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR project with cell type labels transferred from RNA
#' @examples
#' \dontrun{sc_obj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "RNA")}
#' \dontrun{proj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "ATAC")}
#' @export
TransferRNALabels <- function(sc_obj, proj, log_flag = FALSE) {
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 10: Transferring cell type labels from RNA data to ATAC data (true multiome)", log_flag)
  #curated_snRNA_seq_cells <- curated_snRNA_seq_cells[curated_snRNA_seq_cells$cells %in% proj$cellNames,]
  #curated_snRNA_seq_cells <- curated_snRNA_seq_cells[order(match(curated_snRNA_seq_cells$cells,proj$cellNames)),]
  #snRNA_seq_cell_votes <- curated_snRNA_seq_cells$voted_type
  #proj <- addCellColData(ArchRProj = proj, data = snRNA_seq_cell_votes, cells = proj$cellNames, name = "predictedGroup", force = TRUE)
  print_SPEEDI("Step 10: Complete", log_flag)
  gc()
  return(proj)
}
