#' Find overlap between true multiome cells (RNA / ATAC)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param proj ArchR project associated with data
#' @param data_modality data modality to return (in Seurat object or ArchR project)
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object or ArchR project with only those cells that pass QC in the other assay
#' @examples
#' \dontrun{sc_obj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "RNA")}
#' \dontrun{proj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "ATAC")}
#' @export
FindMultiomeOverlap <- function(sc_obj, proj, data_modality = "RNA", output_dir = getwd(), exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  multiome_obj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 9: Finding overlap between RNA and ATAC for multiome data", log_flag)
      print_SPEEDI(paste0("data_modality is: ", data_modality), log_flag)
      if(data_modality == "ATAC") {
        print_SPEEDI("Saving ArchR project (True Multiome)", log_flag = log_flag)
        proj <- ArchR::subsetArchRProject(ArchRProj = proj,
                                          cells = colnames(sc_obj),
                                          force = TRUE,
                                          outputDirectory = "ArchRMultiomeOverlap")
        num_cells <- length(proj$cellNames)
        num_samples <- length(unique(proj$Sample))
        sample_text <- ""
        if(num_samples == 1) {
          sample_text <- paste0("(1 Sample, ", num_cells, " Cells)")
        } else {
          sample_text <- paste0("(", num_samples, " Samples, ", num_cells, " Cells)")
        }
        p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
          ggplot2::ggtitle(paste0("Multiome Overlap for ATAC Data (By Sample) \n ", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18))
        ggplot2::ggsave(filename = paste0(output_dir, "Final_Multiome_Overlap_ATAC_UMAP_by_Sample.png"), plot = p1, device = "png", width = 8, height = 8, units = "in")
      } else if(data_modality == "RNA") {
        sc_obj <- sc_obj[,ArchR::getCellNames(proj)]
        sample_count <- length(unique(sc_obj$sample))
        cell_count <- length(sc_obj$cell_name)
        current_title <- paste0("Multiome Overlap for RNA Data (By Sample) \n (", sample_count, " Samples, ", cell_count, " Cells)")
        print_UMAP_RNA(sc_obj, file_name = "Final_Multiome_Overlap_RNA_UMAP_by_Sample.png",
                       group_by_category = "sample", output_dir = output_dir, title = current_title,
                       log_flag = log_flag)
      } else {
        print_SPEEDI("Error: you have selected an invalid data modality for multiome overlap (you must choose RNA or ATAC).", log_flag)
        exit_code <- 30
        stop()
      }
      print_SPEEDI("Step 9: Complete", log_flag)
      if(data_modality == "ATAC") {
        return(proj)
      } else {
        return(sc_obj)
      }
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running FindMultiomeOverlap() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 27
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(multiome_obj)
}

#' Transfer cell type labels from RNA to ATAC (only works for true multiome since same cells in both assays)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param proj ArchR project associated with data
#' @param RNA_output_dir Path to directory where RNA output will be saved. Defaults to working directory ([getwd()]).
#' @param ATAC_output_dir Path to directory where ATAC output will be saved. Defaults to working directory ([getwd()]).
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR project with cell type labels transferred from RNA
#' @examples
#' \dontrun{sc_obj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "RNA")}
#' \dontrun{proj <- FindMultiomeOverlap(sc_obj, proj, data_modality = "ATAC")}
#' @export
TransferRNALabels <- function(sc_obj, proj, RNA_output_dir = getwd(), ATAC_output_dir = getwd(), exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  proj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 10: Transferring cell type labels from RNA data to ATAC data (true multiome)", log_flag)
      proj <- ArchR::addCellColData(ArchRProj = proj, data = as.character(sc_obj$predicted_celltype_majority_vote), cells = proj$cellNames, name = "Cell_type_voting_RNA", force = TRUE)
      # ATAC
      num_cells <- length(proj$cellNames)
      num_samples <- length(unique(proj$Sample))
      sample_text <- ""
      if(num_samples == 1) {
        sample_text <- paste0("(1 Sample, ", num_cells, " Cells)")
      } else {
        sample_text <- paste0("(", num_samples, " Samples, ", num_cells, " Cells)")
      }
      p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting_RNA", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("Multiome Overlap for ATAC Data (By RNA Majority Vote Cell Type) \n ", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18))
      ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_Multiome_Overlap_ATAC_UMAP_by_RNA_Majority_Vote_Cell_Type_Labels.png"), plot = p1, device = "png", width = 8, height = 8, units = "in")
      p2 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("Multiome Overlap for ATAC Data (By ATAC Majority Vote Cell Type) \n ", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18))
      ggplot2::ggsave(filename = paste0(ATAC_output_dir, "Final_Multiome_Overlap_ATAC_UMAP_by_Original_ATAC_Majority_Vote_Cell_Type_Labels.png"), plot = p2, device = "png", width = 8, height = 8, units = "in")
      ArchR::plotPDF(p1,p2, name = "UMAP_multiome_ATAC_with_RNA_labels_and_original_labels_plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
      # RNA
      sample_count <- length(unique(sc_obj$sample))
      cell_count <- length(sc_obj$cell_name)
      current_title <- paste0("Multiome Overlap for RNA Data (By RNA Majority Vote Cell Type) \n (", sample_count, " Samples, ", cell_count, " Cells)")
      print_UMAP_RNA(sc_obj, file_name = "Final_Multiome_Overlap_RNA_UMAP_by_RNA_Majority_Vote_Cell_Type_Labels.png",
                     group_by_category = "predicted_celltype_majority_vote", output_dir = RNA_output_dir, title = current_title,
                     log_flag = log_flag)
      print_SPEEDI("Step 10: Complete", log_flag)
      return(proj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running TransferRNALabels() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 28
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(proj)
}
