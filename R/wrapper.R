#' Wrapper function for running the SPEEDI pipeline.
#'
#' @description
#' `run_SPEEDI()` performs all steps of the SPEEDI pipeline, including data processing, batch detection and integration, and
#' cell type mapping from your selected reference. Parameters chosen by SPEEDI are recorded in a log file written to your
#' selected output directory (`output_dir`). Relevant individual functions called include:
#'
#' * [Read_h5()]: Reads in input data
#' * [FilterRawData()]: Filters raw data using automatically selected QC thresholds
#' * [InitialProcessing()]: Process filtered data and prepare for batch inferring and integration
#' * [InferBatches()]: Infer batches in data
#' * [IntegrateByBatch()]: Integrate batches together
#' * [VisualizeIntegration()]: Visualize integration and prepare data for marker detection
#' * [LoadReferenceSPEEDI()]: Load reference for mapping onto query data
#' * [MapCellTypes()]: Map cell types from reference onto query data
#' @param tissue Tissue type of input data (used to map cell types via reference). For human data, possible choices include:
#' * `"adipose"`
#' * `"bone marrow"`
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
#'
#' For mouse data, possible choices include:
#' * `"cortex"`
#' * `"custom"` (`reference_file_name` must be provided)
#' @param data_type Type of data being processed. Possible choices include:
#' * `"RNA"`
#' * `"ATAC"` (`tissue` must be set to `"pbmc_full"` or `"custom"`)
#' * `"sample_paired"` (`tissue` must be set to `"pbmc_full"` or `"custom"`)
#' * `"true_multiome"`
#' @param data_path Path to directory where input data are located. Defaults to working directory ([getwd()]).
#' @param reference_dir Path to directory where reference is either already located (see `reference_file_name`) or will be downloaded by [SPEEDI::LoadReferenceSPEEDI()] if necessary. Defaults to working directory ([getwd()]). Note that Azimuth references from [SeuratData] do not require use of `reference_dir`.
#' @param reference_file_name Base name of custom reference file. Should be located inside `reference_dir` and `tissue` should be set to `"custom"`.
#' @param reference_cell_type_attribute If using a Seurat reference object, this parameter captures where the cell type information is stored
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param sample_id_list Vector of sample names (optional - if not provided, will select all samples found recursively in `data_path`).
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param record_doublets Boolean flag to indicate whether we will record doublets in the data (using the [scDblFinder] package). Possible choices are `TRUE` or `FALSE`.
#' @return A Seurat object that has been processed through the SPEEDI pipeline
#' @examples
#' sc_obj <- run_SPEEDI(tissue = "PBMC", species = "human")
#' sc_obj <- run_SPEEDI(tissue = "adipose", data_path = "~/input_data/", reference_dir = "~/reference_dir/", output_dir = "~/adipose_output", species = "human", record_doublets = TRUE)
#' @export
run_SPEEDI <- function(tissue, data_type = "RNA", data_path = getwd(), reference_dir = getwd(), reference_file_name = NULL, reference_cell_type_attribute = "celltype.l2", output_dir = getwd(), sample_id_list = NULL, species = "human", record_doublets = FALSE) {
  # ArchR likes to write some files to the working directory, so we'll set our working directory to output_dir
  # and then reset it to the old working directory once we're done running SPEEDI
  old_wd <- getwd()
  setwd(output_dir)
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
  log_file_name <- paste0(output_dir, log_file_name)
  log_file <- logr::log_open(log_file_name, logdir = FALSE)
  # Error checking
  if((data_type == "ATAC" | data_type == "sample_paired") & (tolower(tissue) != "pbmc_full" & tolower(tissue) != "custom")) {
    print_SPEEDI("Error: You cannot use an Azimuth reference if you are processing ATAC or sample-paired data.", TRUE)
    stop()
  }
  # If there are RNA data, we read those in using Read_RNA, and if there are ATAC data, we read those in using Read_ATAC
  if(data_type != "ATAC") {
    # Read in RNA data, filter data, perform initial processing, infer batches, and integrate by batch
    all_sc_exp_matrices <- Read_RNA(data_path = data_path, sample_id_list = sample_id_list, log_flag = TRUE)
    sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = species,
                                record_doublets = record_doublets, log_file_path = log_file_name, log_flag = TRUE)
    rm(all_sc_exp_matrices)
    sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, log_flag = TRUE)
    sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = TRUE)
    sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
    # Create UMAP of integration (and prep for FindMarkers)
    sc_obj <- VisualizeIntegration(sc_obj = sc_obj, log_flag = TRUE)
  }
  if(data_type != "RNA") {
    # Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
    atac_proj <- Read_ATAC(data_path = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
    atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = TRUE)
    atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = TRUE)
    atac_proj <- IntegrateByBatch_ATAC(proj = atac_proj, log_flag = TRUE)
  }
  # Load reference
  reference <- LoadReferenceSPEEDI(tissue = tissue, species = species, reference_dir = reference_dir,
                                   reference_file_name = reference_file_name, log_flag = TRUE)
  # Map cell types from reference onto query data
  if(data_type != "ATAC") {
    sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                               reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
  }
  if(data_type != "RNA") {
    atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference,
                                   reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
  }
  # Write Seurat object to output directory
  save(sc_obj, file = paste0(log_file_name, ".rds"))
  # Save ArchR project
  saveArchRProject(ArchRProj = atac_proj, load = FALSE)
  setwd(old_wd)
  if(data_type == "RNA") {
    return(sc_obj)
  } else if(data_type == "ATAC") {
    return(atac_proj)
  } else {
    return(list(sc_obj, atac_proj))
  }
}
