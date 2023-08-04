#' Wrapper function for running the SPEEDI pipeline.
#'
#' @description
#' `run_SPEEDI()` performs all steps of the SPEEDI pipeline, including data processing, batch detection and integration, and
#' cell type mapping from your selected reference. Parameters chosen by SPEEDI are recorded in a log file written to your
#' selected output directory (`output_dir`). Relevant individual functions called include:
#'
#' * [Read_RNA()]: Reads in input data (RNA)
#' * [Read_ATAC()]: Reads in input data (ATAC)
#' * [FilterRawData_RNA()]: Filters raw data using automatically selected QC thresholds (RNA)
#' * [FilterRawData_ATAC()]: Filters raw data using automatically selected QC thresholds (ATAC)
#' * [InitialProcessing_RNA()]: Process filtered data and prepare for batch inferring and integration (RNA)
#' * [InitialProcessing_ATAC()]: Process filtered data and prepare for batch inferring and integration (ATAC)
#' * [InferBatches()]: Infer batches in data (both RNA and ATAC)
#' * [IntegrateByBatch_RNA()]: Integrate batches together (RNA)
#' * [IntegrateByBatch_ATAC()]: Integrate batches together (ATAC)
#' * [VisualizeIntegration()]: Visualize integration and prepare data for marker detection (RNA)
#' * [LoadReferenceSPEEDI()]: Load reference for mapping onto query data (both RNA and ATAC)
#' * [MapCellTypes_RNA()]: Map cell types from reference onto query data (RNA)
#' * [MapCellTypes_ATAC()]: Map cell types from reference onto query data (ATAC)
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
#' @return A Seurat object that has been processed through the SPEEDI pipeline
#' @examples
#' \dontrun{sc_obj <- run_SPEEDI(reference_tissue = "PBMC", species = "human")}
#' \dontrun{sc_obj <- run_SPEEDI(reference_tissue = "adipose", data_path = "~/input_data/",
#' reference_dir = "~/reference_dir/", output_dir = "~/adipose_output",
#' species = "human", record_doublets = TRUE)}
#' @export
#' @import ArchR
run_SPEEDI <- function(reference_tissue, data_type = "RNA", species = "human", data_path = getwd(), reference_dir = getwd(), output_dir = getwd(), metadata_df = NULL, reference_file_name = NULL, reference_cell_type_attribute = "celltype.l2", analysis_name = NULL, sample_id_list = NULL, sample_file_paths = NULL, record_doublets = FALSE) {
  SPEEDI_variables <- initialize_SPEEDI(reference_tissue, data_type, species, data_path, reference_dir, output_dir, metadata_df, reference_file_name, reference_cell_type_attribute, analysis_name, sample_id_list, sample_file_paths, record_doublets)
  log_flag <- TRUE
  print_SPEEDI("Beginning SPEEDI Run!", log_flag = log_flag)
  # Load reference
  reference <- LoadReferenceSPEEDI(reference_tissue = SPEEDI_variables$reference_tissue, species = SPEEDI_variables$species, reference_dir = SPEEDI_variables$reference_dir,
                                   reference_file_name = SPEEDI_variables$reference_file_name, log_flag = log_flag)
  if(data_type != "ATAC") {
    # Read in RNA data, filter data, perform initial processing, infer batches, integrate by batch, and process UMAP of integration
    all_sc_exp_matrices <- Read_RNA(data_path = SPEEDI_variables$data_path, sample_id_list = SPEEDI_variables$sample_id_list,
                                    sample_file_paths = SPEEDI_variables$sample_file_paths, log_flag = log_flag)
    sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = SPEEDI_variables$species,
                                record_doublets = SPEEDI_variables$record_doublets, output_dir = SPEEDI_variables$RNA_output_dir,
                                log_file_path = SPEEDI_variables$log_file_name, log_flag = log_flag)
    rm(all_sc_exp_matrices)
    sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = SPEEDI_variables$species, output_dir = SPEEDI_variables$RNA_output_dir,
                                    metadata_df = SPEEDI_variables$metadata_df, log_flag = log_flag)
    sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = log_flag)
    sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = log_flag)
    sc_obj <- VisualizeIntegration(sc_obj = sc_obj, log_flag = log_flag)
  }
  if(data_type != "RNA") {
    # Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
    atac_proj <- Read_ATAC(data_path = SPEEDI_variables$data_path, sample_id_list = SPEEDI_variables$sample_id_list,
                           sample_file_paths = SPEEDI_variables$sample_file_paths, species = SPEEDI_variables$species,
                           log_flag = log_flag)
    atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = log_flag)
    atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = log_flag)
    atac_proj <- IntegrateByBatch_ATAC(proj = atac_proj, log_flag = log_flag)
  }
  # Map cell types from reference onto query data
  if(data_type != "ATAC") {
    sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                               reference_cell_type_attribute = SPEEDI_variables$reference_cell_type_attribute,
                               output_dir = SPEEDI_variables$RNA_output_dir, log_flag = log_flag)
    # If the user provided metadata, we can perform downstream analyses (differential expression, functional module discovery)
    if(!is.null(metadata_df)) {
      run_downstream_analyses_RNA(sc_obj = sc_obj, reference_tissue = SPEEDI_variables$reference_tissue, species = SPEEDI_variables$species, metadata_df = SPEEDI_variables$metadata_df, output_dir = SPEEDI_variables$RNA_output_dir, log_flag = log_flag)
    }
  }
  if(data_type != "RNA") {
    atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference,
                                   reference_cell_type_attribute = SPEEDI_variables$reference_cell_type_attribute, log_flag = log_flag)
  }
  # Write Seurat object to output directory
  if(data_type != "ATAC") {
    print_SPEEDI("Saving Seurat object (RNA)", log_flag = log_flag)
    save(sc_obj, file = paste0(SPEEDI_variables$RNA_output_dir, SPEEDI_variables$analysis_name, ".RNA.rds"))
  }
  # Save ArchR project
  if(data_type != "RNA") {
    print_SPEEDI("Saving ArchR project (ATAC)", log_flag = log_flag)
    ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
  }
  if(data_type == "true_multiome") {
    sc_obj <- FindMultiomeOverlap(sc_obj = sc_obj, proj = atac_proj, data_modality = "RNA", log_flag = log_flag)
    print_SPEEDI("Saving Seurat object (True Multiome)", log_flag = log_flag)
    save(sc_obj, file = paste0(SPEEDI_variables$RNA_output_dir, SPEEDI_variables$analysis_name, ".RNA.multiome.rds"))
    atac_proj <- FindMultiomeOverlap(sc_obj = sc_obj, proj = atac_proj, data_modality = "ATAC", log_flag = log_flag)
    ATAC_multiome_output_dir <- paste0(SPEEDI_variables$ATAC_output_dir, "ArchRMultiomeOutput", "/")
    atac_proj <- TransferRNALabels(sc_obj = sc_obj, proj = atac_proj, log_flag = log_flag)
    print_SPEEDI("Saving ArchR project (True Multiome)", log_flag = log_flag)
    saveArchRProject(ArchRProj = atac_proj, load = FALSE, outputDirectory = ATAC_multiome_output_dir)
  }
  # If any ATAC plots exist, copy them to the ATAC base directory so they're easier for the user to find
  if(data_type != "RNA") {
    atac_plot_files <- list.files(path = SPEEDI_variables$ATAC_output_dir, pattern = "_plots\\.pdf$", recursive = TRUE, full.names = TRUE)
    file.copy(from = atac_plot_files, to = SPEEDI_variables$ATAC_output_dir)
  }
  # Delete Rplots.pdf file if it exists (junk file created by R batch mode)
  if(file.exists(paste0(SPEEDI_variables$output_dir, "Rplots.pdf"))) {
    file.remove(paste0(SPEEDI_variables$output_dir, "Rplots.pdf"))
  }
  setwd(SPEEDI_variables$old_wd)
  print_SPEEDI("SPEEDI Run Complete!", log_flag = log_flag)
  if(data_type == "RNA") {
    return(sc_obj)
  } else if(data_type == "ATAC") {
    return(atac_proj)
  } else {
    return(list(sc_obj, atac_proj))
  }
}
