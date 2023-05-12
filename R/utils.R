#' Scale a vector to range(0,1)
#'
#' @param x Numeric vector
#' @return A scaled vector ranging from 0 to 1
scale_zero_one <- function(x) {(x - min(x))/(max(x) - min(x))}

#' Append a list to a list-of-lists
#'
#' @param lst List
#' @param ... Additional lists
#' @return A list of lists
lappend <- function (lst, ...){ c(lst, list(...))}

#' Parse namespace and function name for do.call
#'
#' @param x Function name (potentially with namespace attached)
#' @return Parsed function
parse_namespace_and_function <- function(x) {
  if(length(grep("::", x)) > 0) {
    parts <- strsplit(x, "::")[[1]]
    getExportedValue(parts[1], parts[2])
  } else {
    x
  }
}

#' Print to console as well as log file (if it's present)
#' @param current_message Message to print
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @param silence_time don't print time in line
#' @return TRUE
print_SPEEDI <- function(current_message, log_flag = FALSE, silence_time = FALSE) {
  if(!silence_time) {
    current_message <- paste0(Sys.time(), ": ", current_message)
  }
  message(current_message)
  if(log_flag) {
    logr::log_print(current_message, console = FALSE, hide_notes = TRUE, blank_after = FALSE)
  }
  return(TRUE)
}

#' Print UMAP for Seurat object (RNA)
#' @param sc_obj Seurat object containing cells for all samples
#' @param file_name File name for saved plot
#' @param group_by_category Category to group samples by
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return TRUE
print_UMAP_RNA <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), log_flag = FALSE) {
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  if(!is.null(group_by_category)) {
    Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
                    label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    Seurat::DimPlot(sc_obj, reduction = "umap", label = TRUE,
                    label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  ggplot2::ggsave(paste0(output_dir, file_name), device = "png", dpi = 300)
  return(TRUE)
}

#' Print heatmap for cell type proportions in Seurat object (RNA)
#' @param sc_obj Seurat object containing cells for all samples
#' @param file_name File name for saved plot
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return TRUE
print_heatmap_cell_type_proportions_RNA <- function(sc_obj, file_name, output_dir = getwd(), log_flag = FALSE) {
  voting_cell_type_proportion <- as.matrix(table(sc_obj$sample, sc_obj$predicted_celltype_majority_vote))
  voting_cell_type_proportion <- apply(voting_cell_type_proportion, 1, function(x){x/sum(x)})
  output.plot <- pheatmap::pheatmap(voting_cell_type_proportion,  cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f", legend = FALSE)
  ggplot2::ggsave(filename = paste0(output_dir, file_name), plot = output.plot, device = "png", width = 8, height = 8, units = "in")
  return(TRUE)
}
