#' Initialize SPEEDI variables to make processing easier
#' While this step isn't strictly necessary for using SPEEDI, it will help make sure
#' that all variables are set properly and will create a log file for SPEEDI
#'
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
#' @return List containing initialized SPEEDI variables
#' @export
initialize_SPEEDI <- function(reference_tissue, data_type = "RNA", species = "human", data_path = getwd(), reference_dir = getwd(), output_dir = getwd(), metadata_df = NULL, reference_file_name = NULL, reference_cell_type_attribute = "celltype.l2", analysis_name = NULL, sample_id_list = NULL, sample_file_paths = NULL, record_doublets = FALSE, exit_with_code = FALSE) {
  if(!is.null(data_path) && !dir.exists(data_path)) {
    print_SPEEDI("Error: Input directory doesn't exist.", log_flag)
    quit_SPEEDI(exit_with_code = exit_with_code, exit_code = 1, log_flag = FALSE)
  }
  if(!is.null(output_dir)) {
    print_SPEEDI("Error: You must provide an output directory (cannot be NULL).", log_flag)
    quit_SPEEDI(exit_with_code = exit_with_code, exit_code = 2, log_flag = FALSE)
  }
  # Normalize paths (in case user provides relative paths)
  if(!is.null(data_path)) {
    data_path <- normalize_dir_path(data_path)
  }
  if(!is.null(reference_dir)) {
    reference_dir <- normalize_dir_path(reference_dir)
  }
  output_dir <- normalize_dir_path(output_dir)
  # Change reference_tissue and species to all lowercase to prevent any issues with casing
  reference_tissue <- tolower(reference_tissue)
  species <- tolower(species)
  # ArchR likes to write some files to the working directory, so we'll set our working directory to output_dir
  # and then reset it to the original working directory once we're done running SPEEDI
  old_wd <- getwd()
  # Create output_dir if it doesn't already exist
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  # Add "/" to end of output_dir if not already present
  last_char_of_output_dir_path <- substr(output_dir, nchar(output_dir), nchar(output_dir))
  if(last_char_of_output_dir_path != "/") {
    output_dir <- paste0(output_dir, "/")
  }
  # Set analysis name
  if(is.null(analysis_name)) {
    analysis_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
    analysis_name <- gsub(":", "-", analysis_name)
  }
  # Update our output dir to be the specific analysis directory
  output_dir <- paste0(output_dir, analysis_name, "/")
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  setwd(output_dir)
  # Create log file
  log_file_info <- create_SPEEDI_log_file(output_dir)
  log_file <- log_file_info[[1]]
  log_file_path <- log_file_info[[2]]
  log_flag <- TRUE
  print_SPEEDI("Log file successfully created", log_flag)
  # If metadata_df is not null, set sample_id_list according to rownames of metadata_df
  if(!is.null(metadata_df)) {
    sample_id_list <- rownames(metadata_df)
  }
  # Output dirs for RNA and ATAC
  RNA_output_dir <- paste0(output_dir, "RNA", "/")
  ATAC_output_dir <- paste0(output_dir, "ATAC", "/")
  if(data_type != "ATAC") {
    if (!dir.exists(RNA_output_dir)) {dir.create(RNA_output_dir)}
  } else if (data_type != "RNA") {
    if (!dir.exists(ATAC_output_dir)) {dir.create(ATAC_output_dir)}
    setwd(ATAC_output_dir)
  }
  # Check user parameters for immediate errors
  preliminary_check_for_SPEEDI_errors(reference_tissue = reference_tissue, data_type = data_type, species = species, data_path = data_path, reference_dir = reference_dir, output_dir = output_dir, metadata_df = metadata_df, reference_file_name = reference_file_name, reference_cell_type_attribute = reference_cell_type_attribute, analysis_name = analysis_name, sample_id_list = sample_id_list, sample_file_paths = sample_file_paths, record_doublets = record_doublets, exit_with_code = exit_with_code, log_flag = log_flag)
  # Return all updated variables in SPEEDI_variables list
  SPEEDI_variables <- list(reference_tissue = reference_tissue, data_type = data_type, species = species, data_path = data_path,
                           reference_dir = reference_dir, output_dir = output_dir, metadata_df = metadata_df,
                           reference_file_name = reference_file_name, reference_cell_type_attribute = reference_cell_type_attribute,
                           analysis_name = analysis_name, sample_id_list = sample_id_list, sample_file_paths = sample_file_paths,
                           record_doublets = record_doublets, RNA_output_dir = RNA_output_dir, ATAC_output_dir = ATAC_output_dir,
                           log_file_path = log_file_path, old_wd = old_wd)
  print_SPEEDI("Updated SPEEDI variables are: ", log_flag)
  for(index in 1:length(SPEEDI_variables)) {
    print_SPEEDI(paste0(names(SPEEDI_variables)[index], ": ", SPEEDI_variables[index]), log_flag)
  }
  return(SPEEDI_variables)
}
