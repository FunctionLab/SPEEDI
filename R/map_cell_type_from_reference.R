#' Load appropriate reference for SPEEDI
#'
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
#' * `"cortex`
#' * `"custom"` (`reference_file_name` must be provided)
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param reference_dir Path to directory where reference is either already located (see `reference_file_name`) or will be downloaded by [SPEEDI::LoadReferenceSPEEDI()] if necessary. Defaults to working directory ([getwd()]). Note that Azimuth references from [SeuratData] do not require use of `reference_dir`.
#' @param reference_file_name Base name of custom reference file. Should be located inside `reference_dir` and `tissue` should be set to `"custom"`.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A reference object
#' @examples
#' reference <- LoadReferenceSPEEDI(tissue = "PBMC", species = "human")
#' reference <- LoadReferenceSPEEDI(tissue = "cortex", species = "mouse")
#' reference <- LoadReferenceSPEEDI(tissue = "custom", species = "human", reference_dir = "~/reference/", reference_file_name = "custom_pbmc_reference.h5")
#' @export
LoadReferenceSPEEDI <- function(tissue, species = "human", reference_dir = getwd(), reference_file_name = NULL, log_flag = FALSE) {
  # Change tissue to all lowercase to prevent any issues with casing
  tissue <- tolower(tissue)
  species <- tolower(species)
  # Add "/" to end of reference path if not already present
  last_char_of_reference_path <- substr(reference_dir, nchar(reference_dir), nchar(reference_dir))
  if(last_char_of_reference_path != "/") {
    reference_dir <- paste0(reference_dir, "/")
  }
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI(paste0("Step 7: loading reference (and installing data if necessary)"), log_flag)
  print_SPEEDI(paste0("tissue is: ", tissue), log_flag)
  print_SPEEDI(paste0("species is: ", species), log_flag)
  print_SPEEDI(paste0("reference_dir (if necessary) is: ", reference_dir), log_flag)
  if(!is.null(reference_file_name)) {
    print_SPEEDI(paste0("reference_file_name is: ", reference_file_name), log_flag)
  }
  if (species == "human") {
    if (tissue == "adipose") {
      SeuratData::InstallData("adiposeref")
      reference <- "adiposeref"
    } else if (tissue == "bone marrow") {
      SeuratData::InstallData("bonemarrowref")
      reference <- "bonemarrowref"
    } else if (tissue == "cortex") {
      SeuratData::InstallData("humancortexref")
      reference <- "humancortexref"
    } else if (tissue == "fetus") {
      SeuratData::InstallData("fetusref")
      reference <- "fetusref"
    } else if (tissue == "heart") {
      SeuratData::InstallData("heartref")
      reference <- "heartref"
    } else if (tissue == "kidney") {
      SeuratData::InstallData("kidneyref")
      reference <- "kidneyref"
    } else if (tissue == "lung") {
      SeuratData::InstallData("lungref")
      reference <- "lungref"
    } else if (tissue == "pancreas") {
      SeuratData::InstallData("pancreasref")
      reference <- "pancreasref"
    } else if (tissue == "pbmc") {
      SeuratData::InstallData("pbmcref")
      reference <- "pbmcref"
    } else if (tissue == "pbmc_full") {
      reference_url <- get_pbmc_reference_url()
      print_SPEEDI(paste0("Downloading PBMC reference from ", reference_url), log_flag)
      # Download PBMC reference if the user doesn't have it
      if(!file.exists(paste0(reference_dir, sub("\\?.*", "", basename(reference_url))))) {
        httr::GET(
          url = reference_url,
          httr::write_disk(paste0(reference_dir, sub("\\?.*", "", basename(reference_url)))),
          httr::verbose()
        ) -> res
      }
      # Load and return PBMC reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, basename(reference_url)))
    } else if (tissue == "tonsil") {
      SeuratData::InstallData("tonsilref")
      reference <- "tonsilref"
    } else if (tissue == "custom") {
      # Load and return reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, reference_file_name))
    } else {
      message(paste0("\nYour tissue ", tissue, " is not valid for the selected species"))
    }
  } else {
    if (tissue == "cortex") {
      SeuratData::InstallData("mousecortexref")
      reference <- "mousecortexref"
    } else if (tissue == "custom") {
      # Load and return reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, reference_file_name))
    } else {
      message(paste0("\nYour tissue ", tissue, " is not valid for the selected species"))
    }
  }
  if(inherits(reference, "character")) {
    print_SPEEDI(paste0("Selected reference based on your tissue is: ", reference), log_flag)
  }
  print_SPEEDI("Step 7: Complete", log_flag)
  return(reference)
}

#' Find mapping anchors between reference and query
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference A Seurat reference object
#' @param data_type string to indicate whether we're analyzing scRNA data or other data
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return Mapping anchors between reference and query
#' @examples
#' anchors <- FindMappingAnchors(sc_obj, reference = custom_reference_seurat_object, data_type = "scRNA")
#' @export
FindMappingAnchors <- function(sc_obj, reference, data_type = "scRNA", log_flag = FALSE) {
  print_SPEEDI("Finding mapping anchors", log_flag)
  # We don't want to recompute residuals if our reference is too different from our data type (e.g., scRNA versus snRNA)
  if(data_type == "scRNA") {
    recompute.residuals.value <- "T"
  } else {
    print_SPEEDI("Not using recompute.residuals for mapping anchors due to non-scRNA query dataset", log_flag)
    recompute.residuals.value <- "F"
  }
  anchors <- Seurat::FindTransferAnchors(reference = reference,
                                           query = sc_obj,
                                           normalization.method = "SCT",
                                           recompute.residuals = recompute.residuals.value,
                                           reference.reduction = "spca")
  print_SPEEDI("Done finding mapping anchors", log_flag)
  return(anchors)
}

#' In each cluster, we vote for majority cell type
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param current_resolution Parameter that indicates resolution for clustering (should not be set too high)
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains majority vote labels
#' @examples
#' sc_obj <- MajorityVote(sc_obj)
#' @export
MajorityVote <- function(sc_obj, current_resolution = 1, log_flag = FALSE) {
  print_SPEEDI("Begin majority voting", log_flag)
  if(Seurat::DefaultAssay(sc_obj) == "integrated") {
    associated_res_attribute <- paste0("integrated_snn_res.", current_resolution)
  } else {
    associated_res_attribute <- paste0("SCT_snn_res.", current_resolution)
  }
  sc_obj <- Seurat::FindNeighbors(sc_obj, reduction = "pca", dims = 1:30)
  sc_obj <- Seurat::FindClusters(sc_obj, resolution = current_resolution)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)

  integrated_snn_res_df <- sc_obj[[associated_res_attribute]]
  integrated_snn_res_cell_names <- rownames(integrated_snn_res_df)
  integrated_snn_res_values <- integrated_snn_res_df[,1]

  cluster.dump <- as.numeric(levels(integrated_snn_res_values))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    print(i)
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(integrated_snn_res_df[cells,]))
    freq.table <- freq.table[order(freq.table$Freq, decreasing = TRUE),]
    freq.table$diff <- abs(c(diff(freq.table$Freq), 0))
    if(nrow(freq.table) > 30) {
      freq.table <- freq.table[1:30,]
    }
    p.values <- outliers::dixon.test(freq.table$diff)$p.value[[1]]
    max.index <- which.max(freq.table$diff)
    clusters <- as.numeric(as.character(freq.table$Var1[1:max.index]))
    levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(clusters)] <- i
    cluster.dump <- cluster.dump[!cluster.dump %in% clusters]
  }

  if (length(cluster.dump) > 0) {
    for (i in cluster.dump) {
      cells <- rownames(subset(integrated_snn_res_df, integrated_snn_res_df[,1] == i,))
      freq.table <- as.data.frame(table(sc_obj$predicted.id[cells]))
      levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(i)] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
    }
  }
  print_SPEEDI("Done with majority voting", log_flag)
  return(sc_obj)
}

#' Choose assay based on whether there are multiple batches (integrated) or only one batch (SCT)
#' @param sc_obj Seurat object containing cells for all samples
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object with default assay appropriately set
#' @examples
#' sc_obj <- SetDefaultAssay(sc_obj)
SetDefaultAssay <- function(sc_obj, log_flag = FALSE) {
  # Assay will be integrated if multiple batches were found - otherwise, we use SCT assay
  if(length(unique(sc_obj$batch)) != 1) {
    Seurat::DefaultAssay(sc_obj) <- "integrated"
  } else {
    Seurat::DefaultAssay(sc_obj) <- "SCT"
  }
  return(sc_obj)
}

#' We always want to use the `predicted.id` column in our Seurat object to determine majority vote.
#' However, if we use [Azimuth::RunAzimuth], there are often predictions made at multiple annotation levels.
#' Depending on the reference, we use `SetPredictedId` to select the proper annotation level for predicted.id.
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference Reference name
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object with predicted.id appropriately set
#' @examples
#' sc_obj <- SetPredictedId(sc_obj, reference = "bonemarrowref")
SetPredictedId <- function(sc_obj, reference, log_flag = FALSE) {
  print_SPEEDI("Choosing appropriate annotation level from reference", log_flag)
  if(reference == "adiposeref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else if(reference == "bonemarrowref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else if (reference == "fetusref") {
    sc_obj$predicted.id <- sc_obj$predicted.annotation.l2
  } else if (reference == "heartref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else if (reference == "humancortexref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else if (reference == "kidneyref") {
    sc_obj$predicted.id <- sc_obj$predicted.annotation.l2
  } else if (reference == "lungref") {
    sc_obj$predicted.id <- sc_obj$predicted.ann_level_3
  } else if (reference == "pancreasref") {
    sc_obj$predicted.id <- sc_obj$predicted.annotation.l1
  } else if (reference == "pbmcref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else if (reference == "tonsilref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l1
  } else if(reference == "mousecortexref") {
    sc_obj$predicted.id <- sc_obj$predicted.celltype.l2
  } else {
    print_SPEEDI("Invalid reference", log_flag)
  }
  return(sc_obj)
}

#' Map cell types for input data
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference Seurat reference object or reference found in [SeuratData]
#' @param data_type String to indicate whether we're analyzing scRNA or snRNA data
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains majority vote labels
#' @examples
#' sc_obj <- MapCellTypes(sc_obj, reference = custom_reference_seurat_object)
#' sc_obj <- MapCellTypes(sc_obj, reference = "adiposeref")
#' @export
MapCellTypes <- function(sc_obj, reference, data_type = "scRNA", log_flag = FALSE) {
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 8: Reference-based cell type mapping", log_flag)
  print_SPEEDI(paste0("reference is: ", reference), log_flag)
  print_SPEEDI(paste0("data_type is: ", data_type), log_flag)
  # Set default assay (to integrated or SCT)
  sc_obj <- SetDefaultAssay(sc_obj)
  # Grab all possible SeuratData references (to make sure that user provided a valid one)
  possible_seuratdata_references <- get_seuratdata_references()
  if(inherits(reference, "Seurat")) {
    anchors <- FindMappingAnchors(sc_obj, reference, data_type, log_flag)
    print_SPEEDI("Mapping reference onto query cells", log_flag)
    sc_obj <- Seurat::MapQuery(anchorset = anchors,
                     query = sc_obj,
                     reference = reference,
                     refdata = "celltype.l2",
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     verbose = TRUE)
    print_SPEEDI("Done mapping reference onto query cells", log_flag)
    sc_obj <- MajorityVote(sc_obj, log_flag = log_flag)
  } else if(inherits(reference, "character") & reference %in% possible_seuratdata_references) {
    print_SPEEDI("Running Azimuth to map reference onto query cells", log_flag)
    sc_obj <- Azimuth::RunAzimuth(query = sc_obj, reference = reference)
    print_SPEEDI("Done running Azimuth to map reference onto query cells", log_flag)
    sc_obj <- SetDefaultAssay(sc_obj)
    sc_obj <- SetPredictedId(sc_obj, reference, log_flag)
    sc_obj <- MajorityVote(sc_obj, log_flag = log_flag)
  } else {
    if(!inherits(reference, "Seurat") & !inherits(reference, "character")) {
      print_SPEEDI(paste0("\nYour reference is not a supported class. It is class ", class(reference), " and should be a Seurat object or a character string."), log_flag)
    } else if(inherits(reference, "character") & !(reference %in% possible_seuratdata_references)) {
      print_SPEEDI(paste0("\nYour reference name (", reference, ") is not valid (it was not found in SeuratData). It should be one of the following: \n"), log_flag)
      print_SPEEDI(paste0(possible_seuratdata_references, collapse = "\n"), log_flag)
    }
  }
  print_SPEEDI("Step 8: Complete", log_flag)
  gc()
  return(sc_obj)
}


