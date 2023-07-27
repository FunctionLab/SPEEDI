#' Load appropriate reference for SPEEDI
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
#' * `"none"` (reference will not be loaded)
#'
#' For mouse data, possible choices include:
#' * `"cortex`
#' * `"custom"` (`reference_file_name` must be provided)
#' #' * `"none"` (reference will not be loaded)
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param reference_dir Path to directory where reference is either already located (see `reference_file_name`) or will be downloaded by [SPEEDI::LoadReferenceSPEEDI()] if necessary. Defaults to working directory ([getwd()]). Note that Azimuth references from [SeuratData] do not require use of `reference_dir`.
#' @param reference_file_name Base name of custom reference file. Should be located inside `reference_dir` and `reference_tissue` should be set to `"custom"`.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A reference object
#' @examples
#' \dontrun{reference <- LoadReferenceSPEEDI(reference_tissue = "PBMC", species = "human")}
#' \dontrun{reference <- LoadReferenceSPEEDI(reference_tissue = "cortex", species = "mouse")}
#' \dontrun{reference <- LoadReferenceSPEEDI(reference_tissue = "custom", species = "human",
#' reference_dir = "~/reference/", reference_file_name = "custom_pbmc_reference.h5")}
#' @export
LoadReferenceSPEEDI <- function(reference_tissue, species = "human", reference_dir = getwd(), reference_file_name = NULL, log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  reference_dir <- normalize_dir_path(reference_dir)
  # Change reference_tissue to all lowercase to prevent any issues with casing
  reference_tissue <- tolower(reference_tissue)
  species <- tolower(species)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI(paste0("Step 1: loading reference (and installing reference data if necessary)"), log_flag)
  print_SPEEDI(paste0("reference_tissue is: ", reference_tissue), log_flag)
  print_SPEEDI(paste0("species is: ", species), log_flag)
  print_SPEEDI(paste0("reference_dir (if necessary) is: ", reference_dir), log_flag)
  if(!is.null(reference_file_name)) {
    print_SPEEDI(paste0("reference_file_name is: ", reference_file_name), log_flag)
  }
  if (species == "human") {
    if (reference_tissue == "adipose") {
      SeuratData::InstallData("adiposeref")
      reference <- "adiposeref"
    } else if (reference_tissue == "bone_marrow") {
      SeuratData::InstallData("bonemarrowref")
      reference <- "bonemarrowref"
    } else if (reference_tissue == "cortex") {
      SeuratData::InstallData("humancortexref")
      reference <- "humancortexref"
    } else if (reference_tissue == "fetus") {
      SeuratData::InstallData("fetusref")
      reference <- "fetusref"
    } else if (reference_tissue == "heart") {
      SeuratData::InstallData("heartref")
      reference <- "heartref"
    } else if (reference_tissue == "kidney") {
      SeuratData::InstallData("kidneyref")
      reference <- "kidneyref"
    } else if (reference_tissue == "lung") {
      SeuratData::InstallData("lungref")
      reference <- "lungref"
    } else if (reference_tissue == "pancreas") {
      SeuratData::InstallData("pancreasref")
      reference <- "pancreasref"
    } else if (reference_tissue == "pbmc") {
      SeuratData::InstallData("pbmcref")
      reference <- "pbmcref"
    } else if (reference_tissue == "pbmc_full") {
      reference_url <- get_pbmc_reference_url()
      # Download PBMC reference if the user doesn't have it
      if(!file.exists(paste0(reference_dir, sub("\\?.*", "", basename(reference_url))))) {
        print_SPEEDI(paste0("Downloading PBMC reference from ", reference_url), log_flag)
        httr::GET(
          url = reference_url,
          httr::write_disk(paste0(reference_dir, sub("\\?.*", "", basename(reference_url)))),
          httr::verbose()
        ) -> res
      }
      # Load and return PBMC reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, basename(reference_url)))
    } else if (reference_tissue == "tonsil") {
      SeuratData::InstallData("tonsilref")
      reference <- "tonsilref"
    } else if (reference_tissue == "custom") {
      # Load and return reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, reference_file_name))
    } else if (reference_tissue == "none") {
      reference <- "none"
    } else {
      message(paste0("\nYour reference tissue ", reference_tissue, " is not valid for the selected species"))
    }
  } else {
    if (reference_tissue == "cortex") {
      SeuratData::InstallData("mousecortexref")
      reference <- "mousecortexref"
    } else if (reference_tissue == "custom") {
      # Load and return reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_dir, reference_file_name))
    } else if (reference_tissue == "none") {
      reference <- "none"
    } else {
      message(paste0("\nYour reference tissue ", reference_tissue, " is not valid for the selected species"))
    }
  }
  if(inherits(reference, "character")) {
    print_SPEEDI(paste0("Selected reference based on reference_tissue is: ", reference), log_flag)
  }
  print_SPEEDI("Step 1: Complete", log_flag)
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
#' \dontrun{anchors <- FindMappingAnchors(sc_obj, reference = custom_reference_seurat_object,
#' data_type = "scRNA")}
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

#' In each cluster, we vote for majority cell type (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param current_resolution Parameter that indicates resolution for clustering (should not be set too high)
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains majority vote labels
#' @examples
#' \dontrun{sc_obj <- MajorityVote_RNA(sc_obj)}
#' @export
MajorityVote_RNA <- function(sc_obj, current_resolution = 2, log_flag = FALSE) {
  print_SPEEDI("Begin majority voting for RNA-seq...", log_flag)
  # TODO: Make this work when we have just one batch
  DefaultAssay(sc_obj) <- "integrated"
  if (is.null(sc_obj@graphs$integrated_snn)) {
    sc_obj <- FindNeighbors(object = sc_obj, reduction = "pca", dims = 1:30)
  } else {
    print_SPEEDI("Neighbors exist. Skipping constructing neighborhood graph...", log_flag)
  }
  sc_obj <- find_clusters_SPEEDI(sc_obj = sc_obj, resolution = current_resolution, log_flag)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)
  votes <- c()
  vote_levels_fix <- as.character(levels(sc_obj$seurat_clusters))
  vote_levels_mod <- as.character(levels(sc_obj$seurat_clusters))

  Idents(sc_obj) <- "seurat_clusters"
  for (i in vote_levels_fix) {
    sub_sc_obj <- subset(sc_obj, idents = i)

    gmeans <- c()
    cell_types <- c()

    for (j in names(prop.table(table(sub_sc_obj$predicted.id))[prop.table(table(sub_sc_obj$predicted.id)) > 0.25])) {
      cell_types <- c(cell_types, j)
      scores <- sub_sc_obj$predicted.id.score[which(sub_sc_obj$predicted.id == j)]
      gmeans <- c(gmeans, exp(mean(log(scores))))
    }

    if(!is.null(cell_types)) {
      vote_levels[vote_levels %in% as.character(i)] <- cell_types[which.max(gmeans)]
    } else {
      vote_levels[vote_levels %in% as.character(i)] <- "Undetermined"
    }
  }

  rare_ct <- names(which(prop.table(table(sc_obj$predicted.id)) < 0.05))
  for (i in vote_levels_fix) {
    sub_sc_obj <- subset(sc_obj, idents = i)
    for (j in rare_ct) {
      if (j %in% sub_sc_obj$predicted.id &
          prop.table(table(sub_sc_obj$predicted.id))[j] > 0.25) {
        vote_levels_mod[which(vote_levels_fix == i)] <- j
      }
    }
  }

  predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(predicted_celltype_majority_vote) <- vote_levels_mod
  predicted_celltype_majority_vote <- as.character(predicted_celltype_majority_vote)
  sc_obj$predicted_celltype_majority_vote <- predicted_celltype_majority_vote

  print_SPEEDI("...End majority voting for RNA-seq", log_flag)
  return(sc_obj)
}

#' In each cluster, we vote for majority cell type (ATAC)
#'
#' @param proj ArchR object containing cells for all samples
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR object which contains majority vote labels
#' @examples
#' \dontrun{proj <- MajorityVote_ATAC(proj)}
#' @export
MajorityVote_ATAC <- function(proj) {
  # TODO: Make this work when we have just one batch
  print_SPEEDI("Begin majority voting for ATAC-seq data...", log_flag)
  tile_reduc <- getReducedDims(ArchRProj = proj, reducedDims = "Harmony", returnMatrix = TRUE)
  tmp <- matrix(rnorm(nrow(tile_reduc) * 3, 10), ncol = nrow(tile_reduc), nrow = 3)
  colnames(tmp) <- rownames(tile_reduc)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
  obj[['Harmony']] <- Seurat::CreateDimReducObject(embeddings=tile_reduc, key='Harmony_', assay='RNA')
  obj <- FindNeighbors(obj, reduction = "Harmony", dims = 1:29)
  obj <- find_clusters_SPEEDI(obj, resolution = 2, log_flag)
  obj <- RunUMAP(obj, reduction = "Harmony", dims = 1:29)
  proj <- addCellColData(
    ArchRProj = proj,
    cells = names(obj$seurat_clusters),
    data = as.character(obj$seurat_clusters),
    name = "seurat_clusters",
    force = TRUE)

  rm(obj)

  Harmony_clusters <- as.factor(proj$seurat_clusters)
  predictedGroup <- proj$predictedGroup
  predictedScore <- proj$predictedScore

  votes <- c()
  vote_levels <- as.character(levels(Harmony_clusters))

  for (i in vote_levels) {
    cluster_cells <- which(proj$seurat_clusters == i)
    sub_predictedGroup <- predictedGroup[cluster_cells]
    sub_predictedScore <- predictedScore[cluster_cells]

    gmeans <- c()
    cell_types <- c()

    for (j in names(prop.table(table(sub_predictedGroup))[prop.table(table(sub_predictedGroup)) > 0.25])) {
      cell_types <- c(cell_types, j)
      cell_type_cells <- which(sub_predictedGroup == j)
      gmeans <- c(gmeans, exp(mean(log(sub_predictedScore[cell_type_cells]))))
    }

    if(!is.null(cell_types)) {
      vote_levels[vote_levels %in% as.character(i)] <- cell_types[which.max(gmeans)]
    } else {
      vote_levels[vote_levels %in% as.character(i)] <- "Undetermined"
    }

  }
  cell_type_voting <- Harmony_clusters
  levels(cell_type_voting) <- vote_levels
  cell_type_voting <- as.character(cell_type_voting)

  proj <- addCellColData(ArchRProj = proj, data = cell_type_voting,
                         cells = proj$cellNames,
                         name = "Cell_type_voting", force = TRUE)

  print_SPEEDI("...End majority voting for ATAC-seq", log_flag)
  return(proj)
}

#' Choose assay based on whether there are multiple batches (integrated) or only one batch (SCT)
#' @param sc_obj Seurat object containing cells for all samples
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object with default assay appropriately set
#' @examples
#' \dontrun{sc_obj <- SetDefaultAssay(sc_obj)}
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
#' \dontrun{sc_obj <- SetPredictedId(sc_obj, reference = "bonemarrowref")}
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

#' Map cell types for input data (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference Seurat reference object or reference found in [SeuratData]
#' @param reference_cell_type_attribute If using a Seurat reference object, this parameter captures where the cell type information is stored
#' @param data_type String to indicate whether we're analyzing scRNA or snRNA data
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains majority vote labels
#' @examples
#' \dontrun{sc_obj <- MapCellTypes_RNA(sc_obj, reference = custom_reference_seurat_object,
#' reference_cell_type_attribute = "Celltype")}
#' \dontrun{sc_obj <- MapCellTypes_RNA(sc_obj, reference = "adiposeref")}
#' @export
MapCellTypes_RNA <- function(sc_obj, reference, reference_cell_type_attribute = "celltype.l2", data_type = "scRNA", output_dir = getwd(), log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  output_dir <- normalizePath(output_dir, "/")
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 8: Reference-based cell type mapping (RNA)", log_flag)
  if(inherits(reference, "character")) {
    print_SPEEDI(paste0("reference is: ", reference), log_flag)
  }
  if(!is.null(reference_cell_type_attribute)) {
    print_SPEEDI(paste0("reference_cell_type_attribute is: ", reference_cell_type_attribute), log_flag)
  }
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
                     refdata = reference_cell_type_attribute,
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     verbose = TRUE)
    print_SPEEDI("Done mapping reference onto query cells", log_flag)
    sc_obj <- MajorityVote_RNA(sc_obj, log_flag = log_flag)
  } else if(inherits(reference, "character") & reference %in% possible_seuratdata_references) {
    print_SPEEDI("Running Azimuth to map reference onto query cells", log_flag)
    sc_obj <- Azimuth::RunAzimuth(query = sc_obj, reference = reference)
    print_SPEEDI("Done running Azimuth to map reference onto query cells", log_flag)
    sc_obj <- SetDefaultAssay(sc_obj)
    sc_obj <- SetPredictedId(sc_obj, reference, log_flag)
    sc_obj <- MajorityVote_RNA(sc_obj, log_flag = log_flag)
  } else if(inherits(reference, "character") && reference == "none") {
    print_SPEEDI("Not performing reference mapping because no reference was provided", log_flag)
  } else {
    if(!inherits(reference, "Seurat") & !inherits(reference, "character")) {
      print_SPEEDI(paste0("\nYour reference is not a supported class. It is class ", class(reference), " and should be a Seurat object or a character string."), log_flag)
    } else if(inherits(reference, "character") & !(reference %in% possible_seuratdata_references)) {
      print_SPEEDI(paste0("\nYour reference name (", reference, ") is not valid (it was not found in SeuratData). It should be one of the following: \n"), log_flag)
      print_SPEEDI(paste0(possible_seuratdata_references, collapse = "\n"), log_flag)
    }
  }
  if(inherits(reference, "Seurat") || (inherits(reference, "character") && reference != "none")) {
    print_SPEEDI("Printing final UMAPs", log_flag)
    print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Majority_Vote_Cell_Type.png",
               group_by_category = "predicted_celltype_majority_vote", output_dir = output_dir,
               log_flag = log_flag)
    print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Cluster.png",
               group_by_category = "seurat_clusters", output_dir = output_dir,
               log_flag = log_flag)
    print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Raw_Predicted_Cell_Type.png",
               group_by_category = "predicted.id", output_dir = output_dir,
               log_flag = log_flag)
    print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Sample.png",
               group_by_category = "sample", output_dir = output_dir,
               log_flag = log_flag)
    print_heatmap_cell_type_proportions_RNA(sc_obj, file_name = "Final_RNA_Cell_Type_Proportion_Heatmap.png",
                                            output_dir = output_dir, log_flag = log_flag)
  }
  print_SPEEDI("Step 8: Complete", log_flag)
  gc()
  return(sc_obj)
}

#' Map cell types for input data (ATAC)
#'
#' @param proj ArchR project containing cells for all samples
#' @param reference Seurat reference object
#' @param reference_cell_type_attribute If using a Seurat reference object, this parameter captures where the cell type information is stored
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR object which contains majority vote labels
#' @examples
#' \dontrun{sc_obj <- MapCellTypes_ATAC(proj, reference = custom_reference_seurat_object)}
#' @export
MapCellTypes_ATAC <- function(proj, reference, reference_cell_type_attribute = "celltype.l2", log_flag = FALSE) {
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 8: Reference-based cell type mapping (ATAC)", log_flag)
  if(inherits(reference, "character")) {
    print_SPEEDI(paste0("reference is: ", reference), log_flag)
  }
  if(!is.null(reference_cell_type_attribute)) {
    print_SPEEDI(paste0("reference_cell_type_attribute is: ", reference_cell_type_attribute), log_flag)
  }
  if(inherits(reference, "character") && reference == "none") {
    print_SPEEDI("Not performing reference mapping because no reference was provided", log_flag)
  } else {
    print_SPEEDI("Adding gene integration matrix into ArchR project using reference", log_flag)
    # If we only had one batch, then we just use IterativeLSI - otherwise, we use Harmony
    if(length(unique(proj$Batch)) == 1) {
      reducedDims_param <- "IterativeLSI"
    } else {
      reducedDims_param <- "Harmony"
    }
    if(Seurat::DefaultAssay(reference) == "SCT") {
      normalization_method <- "SCT"
    } else {
      normalization_method <- "LogNormalize"
    }
    proj <- addGeneIntegrationMatrix_SPEEDI(
      ArchRProj = proj,
      useMatrix = "GeneScoreMatrix",
      matrixName = "GeneIntegrationMatrix",
      reducedDims = reducedDims_param,
      seRNA = reference,
      dimsToUse = 2:30,
      addToArrow = FALSE,
      groupRNA = reference_cell_type_attribute,
      nameCell = "predictedCell",
      nameGroup = "predictedGroup",
      nameScore = "predictedScore",
      normalization.method = normalization_method,
      force = TRUE
    )
    print_SPEEDI("Done adding gene integration matrix into ArchR project using reference", log_flag)
    pal <- paletteDiscrete(values = proj$predictedGroup)
    p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
    ArchR::plotPDF(p1, name = "UMAP_after_Initial_Cell_Type_Reference_Mapping_plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
    print_SPEEDI("Performing majority voting", log_flag)
    # We have to perform majority voting with a different cluster attribute if Harmony was not run
    # (due to only having one batch)
    cluster_variable_name <- NULL
    if(reducedDims_param == "Harmony") {
      cluster_variable_name <- "Harmony_clusters"
      cM <- as.matrix(ArchR::confusionMatrix(proj$Harmony_clusters, proj$predictedGroup))
      Cell_type_voting <- proj$Harmony_clusters
      pre_cluster <- rownames(cM)
      max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]
      for (m in c(1:length(pre_cluster))){
        idxSample <- which(proj$Harmony_clusters == pre_cluster[m])
        Cell_type_voting[idxSample] <- max_celltype[m]
      }
    } else {
      cluster_variable_name <- "Clusters"
      cM <- as.matrix(ArchR::confusionMatrix(proj$Clusters, proj$predictedGroup))
      Cell_type_voting <- proj$Clusters
      pre_cluster <- rownames(cM)
      max_celltype <- colnames(cM)[apply(cM, 1 , which.max)]
      for (m in c(1:length(pre_cluster))){
        idxSample <- which(proj$Clusters == pre_cluster[m])
        Cell_type_voting[idxSample] <- max_celltype[m]
      }
    }
    proj <- ArchR::addCellColData(ArchRProj = proj, data = Cell_type_voting, cells = proj$cellNames, name = "Cell_type_voting", force = TRUE)
    print_SPEEDI("Done performing majority voting", log_flag)
    pal <- paletteDiscrete(values = proj$Cell_type_voting)
    p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Cell_type_voting", embedding = "UMAP", pal = pal, force = TRUE, keepAxis = TRUE)
    p2 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = cluster_variable_name, embedding = "UMAP", force = TRUE, keepAxis = TRUE)
    p3 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
    p4 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE)
    ArchR::plotPDF(p1,p2,p3,p4, name = "UMAP_after_Final_Cell_Type_Majority_Voting_plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
  }
  print_SPEEDI("Step 8: Complete", log_flag)
  gc()
  return(proj)
}

