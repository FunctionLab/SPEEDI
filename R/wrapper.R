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
#' @param record_doublets Boolean flag to indicate whether we will record doublets in the data (using the [scDblFinder] package). Possible choices are `TRUE` or `FALSE`.
#' @return A Seurat object that has been processed through the SPEEDI pipeline
#' @examples
#' \dontrun{sc_obj <- run_SPEEDI(reference_tissue = "PBMC", species = "human")}
#' \dontrun{sc_obj <- run_SPEEDI(reference_tissue = "adipose", data_path = "~/input_data/",
#' reference_dir = "~/reference_dir/", output_dir = "~/adipose_output",
#' species = "human", record_doublets = TRUE)}
#' @export
#' @import ArchR
run_SPEEDI <- function(reference_tissue, data_type = "RNA", species = "human", data_path = getwd(), reference_dir = getwd(), output_dir = getwd(), metadata_df = NULL, reference_file_name = NULL, reference_cell_type_attribute = "celltype.l2", analysis_name = NULL, sample_id_list = NULL, record_doublets = FALSE) {
  # Normalize paths (in case user provides relative paths)
  data_path <- normalize_dir_path(data_path)
  reference_dir <- normalize_dir_path(reference_dir)
  output_dir <- normalize_dir_path(output_dir)
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
  log_file_name <- paste0(gsub(" ", "_", Sys.time()), "_SPEEDI")
  log_file_name <- gsub(":", "-", log_file_name)
  log_file_name <- paste0(output_dir, log_file_name)
  log_file <- logr::log_open(log_file_name, logdir = FALSE)
  print_SPEEDI("Beginning SPEEDI Run!", log_flag = TRUE)
  # Error checking
  if((data_type == "ATAC" | data_type == "sample_paired") & (tolower(reference_tissue) != "pbmc_full"
                                                             & tolower(reference_tissue) != "custom") & tolower(reference_tissue) != "none") {
    print_SPEEDI("Error: You cannot use an Azimuth reference if you are processing ATAC or sample-paired data.", TRUE)
    stop()
  }
  # Load reference
  reference <- LoadReferenceSPEEDI(reference_tissue = reference_tissue, species = species, reference_dir = reference_dir,
                                   reference_file_name = reference_file_name, log_flag = TRUE)
  # If metadata_df is not null, set sample_id_list according to rownames of metadata_df
  if(!is.null(metadata_df)) {
    sample_id_list <- rownames(metadata_df)
  }
  # Output dirs for RNA and ATAC
  RNA_output_dir <- paste0(output_dir, "RNA", "/")
  ATAC_output_dir <- paste0(output_dir, "ATAC", "/")

  if(data_type != "ATAC") {
    if (!dir.exists(RNA_output_dir)) {dir.create(RNA_output_dir)}
    # Read in RNA data, filter data, perform initial processing, infer batches, integrate by batch, and process UMAP of integration
    all_sc_exp_matrices <- Read_RNA(data_path = data_path, sample_id_list = sample_id_list, log_flag = TRUE)
    sc_obj <- FilterRawData_RNA(all_sc_exp_matrices = all_sc_exp_matrices, species = species,
                                record_doublets = record_doublets, output_dir = RNA_output_dir,
                                log_file_path = log_file_name, log_flag = TRUE)
    rm(all_sc_exp_matrices)
    sc_obj <- InitialProcessing_RNA(sc_obj = sc_obj, species = species, metadata_df = metadata_df, log_flag = TRUE)
    sc_obj <- InferBatches(sc_obj = sc_obj, log_flag = TRUE)
    sc_obj <- IntegrateByBatch_RNA(sc_obj = sc_obj, log_flag = TRUE)
    sc_obj <- VisualizeIntegration(sc_obj = sc_obj, log_flag = TRUE)
  }

  if(data_type != "RNA") {
    if (!dir.exists(ATAC_output_dir)) {dir.create(ATAC_output_dir)}
    setwd(ATAC_output_dir)
    # Read in ATAC data, filter data, perform initial processing, infer batches, and integrate by batch
    atac_proj <- Read_ATAC(data_path = data_path, sample_id_list = sample_id_list, species = species, log_flag = TRUE)
    atac_proj <- FilterRawData_ATAC(proj = atac_proj, log_flag = TRUE)
    atac_proj <- InitialProcessing_ATAC(proj = atac_proj, log_flag = TRUE)
    atac_proj <- IntegrateByBatch_ATAC(proj = atac_proj, log_flag = TRUE)
  }

  # Map cell types from reference onto query data
  if(data_type != "ATAC") {
    sc_obj <- MapCellTypes_RNA(sc_obj = sc_obj, reference = reference,
                               reference_cell_type_attribute = reference_cell_type_attribute,
                               output_dir = RNA_output_dir, log_flag = TRUE)
    # If the user provided metadata, we can perform downstream analyses (differential expression, functional module discovery)
    if(!is.null(metadata_df)) {
      # We run differential expression on each metadata attribute provided by the user
      differential_expression_results <- RunDE_RNA(sc_obj, metadata_df, output_dir = RNA_output_dir, log_flag = TRUE)
      # Next, if species is human, we will run functional module discovery using the DE results
      # TODO: Pick appropriate networks depending on cell type
      if(species == "human") {
        # TODO: Provide this list to user in return statement?
        FMD_results <- list()
        index <- 1
        for(current_de_result in differential_expression_results) {
          for(current_cell_type in unique(current_de_result$Cell_Type)) {
            # Subset to DE results for our current cell type
            cell_specific_de_result <- current_de_result[current_de_result$Cell_Type == current_cell_type,]
            # Grab HB networks associated with cell type and reference tissue
            hb_networks <- grab_hb_networks(current_cell_type, reference_tissue)
            # We grab high FC genes (positive FC) - if there are at least 20 genes, we can do FMD
            high_genes <- cell_specific_de_result[cell_specific_de_result$sc_log2FC > 0.1,]$Gene_Name
            # We also grab highly negative FC genes (opposite direction) - if there are at least 20 genes, we can do FMD
            low_genes <- cell_specific_de_result[cell_specific_de_result$sc_log2FC < -0.1,]$Gene_Name
            # Run FMD for each network
            for(network in hb_networks) {
              # Run FMD for high genes
              FMD_result_high <- run_fmd_wrapper(high_genes, network, RNA_output_dir, current_cell_type, unique(current_de_result$metadata_attribute), "high", log_flag = TRUE)
              if(!is.null(FMD_result_high)) {
                 FMD_results[[index]] <- FMD_result_high
                 index <- index + 1
              }
              # Run FMD for low genes
              FMD_result_low <- run_fmd_wrapper(low_genes, network, RNA_output_dir, current_cell_type, unique(current_de_result$metadata_attribute), "low", log_flag = TRUE)
              if(!is.null(FMD_result_low)) {
                FMD_results[[index]] <- FMD_result_low
                index <- index + 1
              }
            }
          }
        }
      }
    }
  }
  if(data_type != "RNA") {
    atac_proj <- MapCellTypes_ATAC(proj = atac_proj, reference = reference,
                                   reference_cell_type_attribute = reference_cell_type_attribute, log_flag = TRUE)
  }
  # Write Seurat object to output directory
  if(data_type != "ATAC") {
    print_SPEEDI("Saving Seurat object (RNA)", log_flag = TRUE)
    save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".RNA.rds"))
  }
  # Save ArchR project
  if(data_type != "RNA") {
    print_SPEEDI("Saving ArchR project (ATAC)", log_flag = TRUE)
    ArchR::saveArchRProject(ArchRProj = atac_proj, load = FALSE)
  }
  if(data_type == "true_multiome") {
    sc_obj <- FindMultiomeOverlap(sc_obj = sc_obj, proj = atac_proj, data_modality = "RNA", log_flag = TRUE)
    print_SPEEDI("Saving Seurat object (True Multiome)", log_flag = TRUE)
    save(sc_obj, file = paste0(RNA_output_dir, analysis_name, ".RNA.multiome.rds"))
    atac_proj <- FindMultiomeOverlap(sc_obj = sc_obj, proj = atac_proj, data_modality = "ATAC", log_flag = TRUE)
    ATAC_multiome_output_dir <- paste0(ATAC_output_dir, "ArchRMultiomeOutput", "/")
    atac_proj <- TransferRNALabels(sc_obj = sc_obj, proj = atac_proj, log_flag = TRUE)
    print_SPEEDI("Saving ArchR project (True Multiome)", log_flag = TRUE)
    saveArchRProject(ArchRProj = atac_proj, load = FALSE, outputDirectory = ATAC_multiome_output_dir)
  }
  # If any ATAC plots exist, copy them to the ATAC base directory so they're easier for the user to find
  if(data_type != "RNA") {
    atac_plot_files <- list.files(path = ATAC_output_dir, pattern = "_plots\\.pdf$", recursive = TRUE, full.names = TRUE)
    file.copy(from = atac_plot_files, to = ATAC_output_dir)
  }
  # Delete Rplots.pdf file if it exists (junk file created by R batch mode)
  if(file.exists(paste0(output_dir, "Rplots.pdf"))) {
    file.remove(paste0(output_dir, "Rplots.pdf"))
  }
  setwd(old_wd)
  print_SPEEDI("SPEEDI Run Complete!", log_flag = TRUE)
  if(data_type == "RNA") {
    return(sc_obj)
  } else if(data_type == "ATAC") {
    return(atac_proj)
  } else {
    return(list(sc_obj, atac_proj))
  }
}
