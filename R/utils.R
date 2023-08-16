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

#' Normalize directory path for processing
#'
#' @param path Directory path
#' @return Normalized directory path
#' @export
normalize_dir_path <- function(path) {
  path <- normalizePath(path, "/")
  # Add "/" to end of path if not already present
  last_char_of_path <- substr(path, nchar(path), nchar(path))
  if(last_char_of_path != "/") {
    path <- paste0(path, "/")
  }
  return(path)
}

#' Creates log file for SPEEDI
#' @param output_dir Output directory
#' @param log_file_name Log file name
#' @return Open log file
#' @export
create_SPEEDI_log_file <- function(output_dir = getwd(), log_file_name = NULL) {
  # Create log file name based on timestamp
  if(is.null(log_file_name)) {
    log_file_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
    log_file_name <- gsub(":", "-", log_file_name)
    log_file_name <- paste0(output_dir, log_file_name)
  }
  log_file <- logr::log_open(log_file_name, logdir = FALSE)
  return(list(log_file, log_file_name))
}

#' Preliminary check for SPEEDI errors
#' @param reference_tissue Reference tissue type (used to map cell types via reference). For human data, possible choices include:
#' * `"adipose"`
#' * `"bone_marrow"`
#' * `"cortex"`
#' * `"fetus"`
#' * `"heart"`
#' * `"kidney"`
#' * `"lung"`
#' * `"pancreas"`
#' * `"pbmc"` (uses [Azimuth] reference)
#' * `"pbmc_full"` (downloads more complete PBMC reference - should provide better mapping quality but also quite large (~8 GB) to load into memory)
#' * `"tonsil"`
#' * `"custom"` (`reference_file_name` must be provided)
#' * `"none"` (cell mapping from reference will not be performed)
#'
#' For mouse data, possible choices include:
#' * `"cortex"`
#' * `"custom"` (`reference_file_name` must be provided)
#' * `"none"` (cell mapping from reference will not be performed)
#' @param data_type Type of data being processed. Possible choices include:
#' * `"RNA"`
#' * `"ATAC"` (`reference_tissue` must be set to `"pbmc_full"`, `"custom"`, or `"none"`)
#' * `"sample_paired"` (reference mapping for ATAC will not be performed if an Azimuth reference is selected)
#' * `"true_multiome"` (reference mapping for full ATAC will not be performed if an Azimuth reference is selected, but RNA cell types will be mapped over to ATAC for overlapping cells that pass QC)
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param data_path Path to directory where input data are located. Defaults to working directory ([getwd()]).
#' @param reference_dir Path to directory where reference is either already located (see `reference_file_name`) or will be downloaded by [SPEEDI::LoadReferenceSPEEDI()] if necessary. Defaults to working directory ([getwd()]). Note that Azimuth references from [SeuratData] do not require use of `reference_dir`.
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param metadata_df Dataframe containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param reference_file_name Base name of custom reference file. Should be located inside `reference_dir` and `reference_tissue` should be set to `"custom"`.
#' @param reference_cell_type_attribute If using a Seurat reference object, this parameter captures where the cell type information is stored
#' @param analysis_name Name used to create subdirectory in `output_dir` for current analysis run. Directory will be created if it doesn't already exist.
#' @param sample_id_list Vector of sample names (optional - if not provided, will select all samples found recursively in `data_path`).
#' @param sample_file_paths Vector of sample file paths (optional - if not provided, will select all samples found recursively in `data_path`). If using Market Exchange (MEX) Format (matrix.mtx / barcodes.tsv / features.tsv or genes.tsv), please provide a full set of sample paths for only one type of file (e.g., `"c("sample1/matrix.mtx", "sample2/matrix.mtx"`"). If this argument is used, `sample_id_list` is required and should be written in the same order as the sample file paths.
#' @param record_doublets Boolean flag to indicate whether we will record doublets in the data (using the [scDblFinder] package). Possible choices are `TRUE` or `FALSE`.
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag Boolean flag to indicate whether we're also printing to log file
#' @return TRUE
#' @export
preliminary_check_for_SPEEDI_errors <- function(reference_tissue, data_type = "RNA", species = "human", data_path = getwd(), reference_dir = getwd(), output_dir = getwd(), metadata_df = NULL, reference_file_name = NULL, reference_cell_type_attribute = "celltype.l2", analysis_name = NULL, sample_id_list = NULL, sample_file_paths = NULL, record_doublets = FALSE, exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  possible_references <- get_references()
  azimuth_references <- get_azimuth_references()
  if(!(tolower(reference_tissue) %in% possible_references)) {
    print_SPEEDI("Error: Reference is not supported. It should be one of the following: \n", log_flag)
    print_SPEEDI(paste0(possible_references, collapse = "\n"), silence_time = TRUE, log_flag = log_flag)
    exit_code <- 4
  }
  if(data_type != "RNA" && data_type != "ATAC" && data_type != "sample_paired" && data_type != "true_multiome") {
    print_SPEEDI("Error: Data type is not supported (must be RNA, ATAC, sample_paired, or true_multiome).", log_flag = log_flag)
    exit_code <- 5
  }
  if(tolower(species) != "human" && tolower(species) != "mouse") {
    print_SPEEDI("Error: Species is not supported (must be human or mouse).", log_flag = log_flag)
    exit_code <- 6
  }
  if(is.null(data_path) & is.null(sample_file_paths)) {
    print_SPEEDI("Error: You must provide a value for \"data_path\" if you do not provide a value for \"sample_file_paths\".", log_flag = log_flag)
    exit_code <- 7
  }
  if(!is.null(reference_file_name) && is.null(reference_dir)) {
    print_SPEEDI("Error: If reference file name is provided, reference dir must also be provided.", log_flag = log_flag)
    exit_code <- 8
  }
  if(!is.null(reference_file_name) && !file.exists(paste0(reference_dir, reference_file_name))) {
    print_SPEEDI("Error: Reference file was not found in reference dir.", log_flag = log_flag)
    exit_code <- 9
  }
  if(!is.null(sample_file_paths) & is.null(sample_id_list)) {
    print_SPEEDI("Error: You must provide a value for \"sample_id_list\" if you provide a value for \"sample_file_paths\".", log_flag = log_flag)
    exit_code <- 10
  }
  if(!inherits(record_doublets, "logical")) {
    print_SPEEDI("Error: \"record_doublets\" must be TRUE or FALSE.", log_flag = log_flag)
    exit_code <- 11
  }
  if((data_type == "ATAC" | data_type == "sample_paired") & (tolower(reference_tissue) %in% azimuth_references)) {
    print_SPEEDI("Error: You cannot use an Azimuth reference if you are processing ATAC or sample-paired data.", log_flag = log_flag)
    exit_code <- 12
  }
  if(data_type != "RNA" && grepl(" ", output_dir, fixed = TRUE)) {
    print_SPEEDI("Error: ArchR (used for ATAC processing) does not allow spaces in your output directory path.", log_flag = log_flag)
    exit_code <- 31
  }
  return(exit_code)
}

#' Print to console as well as log file (if it's present)
#' @param current_message Message to print
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @param silence_time don't print time in line
#' @return TRUE
#' @export
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

#' Quits SPEEDI (either with exit code if exit_with_code == TRUE or with stop())
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param exit_code exit code
#' @param log_flag Boolean to indicate whether we're also printing to log file
#' @return TRUE
#' @export
quit_SPEEDI <- function(exit_with_code, exit_code = 0, log_flag = FALSE) {
  logr::log_close()
  if(exit_with_code) {
    quit(status = exit_code)
  } else {
    stop_quietly()
  }
  return(TRUE) # Never reached
}

#' Stops SPEEDI quietly (so we don't write superfluous error messages)
#' @return TRUE
#' @export
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
  return(TRUE)
}

#' Print UMAP for Seurat object (RNA)
#' @param sc_obj Seurat object containing cells for all samples
#' @param file_name File name for saved plot
#' @param group_by_category Category to group samples by
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param title Title of plot (if NULL, automatically generated)
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return TRUE
print_UMAP_RNA <- function(sc_obj, file_name, group_by_category = NULL, output_dir = getwd(), title = NULL, log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  current_title <- ""
  if(is.null(title)) {
    current_title <- paste0("RNA Data Integration \n (", sample_count, " Samples, ", cell_count, " Cells)")
  } else {
    current_title <- title
  }
  if(!is.null(group_by_category)) {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", group.by = group_by_category, label = TRUE,
                    label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  } else {
    p <- Seurat::DimPlot(sc_obj, reduction = "umap", label = TRUE,
                    label.size = 3, repel = TRUE, raster = FALSE) +
      ggplot2::labs(title = current_title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  ggplot2::ggsave(paste0(output_dir, file_name), plot = p, device = "png", dpi = 300)
  return(TRUE)
}

#' Print heatmap for cell type proportions in Seurat object (RNA)
#' @param sc_obj Seurat object containing cells for all samples
#' @param file_name File name for saved plot
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return TRUE
print_heatmap_cell_type_proportions_RNA <- function(sc_obj, file_name, output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalize_dir_path(output_dir)
  voting_cell_type_proportion <- as.matrix(table(sc_obj$sample, sc_obj$predicted_celltype_majority_vote))
  voting_cell_type_proportion <- apply(voting_cell_type_proportion, 1, function(x){x/sum(x)})
  sample_count <- length(unique(sc_obj$sample))
  cell_count <- length(sc_obj$cell_name)
  output.plot <- pheatmap::pheatmap(voting_cell_type_proportion,  cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.3f", legend = FALSE) +
    ggplot2::ggtitle(paste0("Cell Type Proportions \n Integrated RNA Data \n (", sample_count, " Samples, ", cell_count, " Cells)")) + ggplot2::theme(plot.title = ggplot2::element_text(size=18))
  ggplot2::ggsave(filename = paste0(output_dir, file_name), plot = output.plot, device = "png", width = 8, height = 8, units = "in")
  return(TRUE)
}

#' First, try Leiden algorithm for clustering. If that doesn't work, then try Louvain
#' @param sc_obj Seurat object containing cells for all samples
#' @param resolution Resolution for clustering
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return Seurat object with clustering complete
#' @export
find_clusters_SPEEDI <- function(sc_obj, resolution, log_flag = FALSE) {
  sc_obj <- tryCatch(
    {
      print_SPEEDI("Trying to use Leiden algorithm for clustering", log_flag)
      Seurat::FindClusters(object = sc_obj, resolution = resolution, algorithm = 4, method='igraph', random.seed = get_speedi_seed())
    },
    error=function(cond) {
      print_SPEEDI("Error using Leiden algorithm for clustering", log_flag)
      print_SPEEDI(cond, log_flag)
      print_SPEEDI("Trying Louvain instead", log_flag)
      return(Seurat::FindClusters(object = sc_obj, resolution = resolution, algorithm = 2, random.seed = get_speedi_seed()))
    },
    finally={
      print_SPEEDI("Clustering complete", log_flag)
    }
  )
  return(sc_obj)
}
