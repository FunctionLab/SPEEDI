#' Load appropriate reference
#'
#' @param tissue tissue associated with samples
#' @param human flag to indicate whether sample is human or mouse
#' @param reference_path path to base directory for reference (or maybe full path to reference?)
#' @return A reference object
#' @export
LoadReference <- function(tissue, human, reference_path = getwd()) {
  # Add "/" to end of reference path if not already present
  last_char_of_reference_path <- substr(reference_path, nchar(reference_path), nchar(reference_path))
  if(last_char_of_reference_path != "/") {
    reference_path <- paste0(reference_path, "/")
  }
  message("Loading reference")
  message(paste0("Installing data for ", tissue, " reference if necessary"))
  if (human) {
    if (tissue == "Adipose") {
      SeuratData::InstallData("adiposeref")
    } else if (tissue == "Bone Marrow") {
      SeuratData::InstallData("bonemarrowref")
    } else if (tissue == "Fetus") {
      SeuratData::InstallData("fetusref")
    } else if (tissue == "Heart") {
      SeuratData::InstallData("heartref")
    } else if (tissue == "Cortex") {
      SeuratData::InstallData("humancortexref")
    } else if (tissue == "Kidney") {
      SeuratData::InstallData("kidneyref")
    } else if (tissue == "Lung") {
      SeuratData::InstallData("lungref")
    } else if (tissue == "Pancreas") {
      SeuratData::InstallData("pancreasref")
    } else if (tissue == "PBMC") {
      reference_url <- get_pbmc_reference_url()
      # Download PBMC reference if the user doesn't have it
      if(!file.exists(paste0(reference_path, basename(reference_url)))) {
        httr::GET(
          url = reference_url,
          httr::write_disk(paste0(reference_path, basename(reference_url))),
          httr::verbose()
        ) -> res
      }
      # Load and return PBMC reference
      reference <- SeuratDisk::LoadH5Seurat(paste0(reference_path, basename(reference_url)))
      return(reference)
    } else if (tissue == "Tonsil") {
      SeuratData::InstallData("tonsilref")
    }
  } else {
    if (tissue == "Cortex") {
      SeuratData::InstallData("mousecortexref")
    }
  }
}

#' Find mapping anchors between reference and query
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference A Seurat reference object
#' @param data_type string to indicate whether we're analyzing scRNA or snRNA data
#' @return Mapping anchors between reference and query
#' @export
FindMappingAnchors <- function(sc_obj, reference, data_type = "scRNA") {
  if(length(unique(sc_obj$batch)) != 1) {
    Seurat::DefaultAssay(sc_obj) <- "integrated"
  }
  # We don't want to recompute residuals if our reference is too different from our data type (e.g., scRNA versus snRNA)
  if(data_type == "scRNA") {
    recompute.residuals.value <- "T"
  } else {
    recompute.residuals.value <- "F"
  }
  anchors <- Seurat::FindTransferAnchors(reference = reference,
                                 query = sc_obj,
                                 normalization.method = "SCT",
                                 recompute.residuals = recompute.residuals.value,
                                 reference.reduction = "spca")
  return(anchors)
}

#' In each cluster, vote for majority cell type
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param current_resolution parameter that indicates current resolution for clustering
#' @return A Seurat object which contains majority vote labels
#' @export
MajorityVote <- function(sc_obj, current_resolution = 1) {
  message("Begin majority voting...")
  if(length(unique(sc_obj$batch)) != 1) {
    Seurat::DefaultAssay(sc_obj) <- "integrated"
    associated_res_attribute <- paste0("integrated_snn_res.", current_resolution)
  } else {
    Seurat::DefaultAssay(sc_obj) <- "SCT"
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


  message("...End majority voting")
  return(sc_obj)
}

#' Map cell types for input data
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param reference Seurat reference object
#' @param data_type string to indicate whether we're analyzing scRNA or snRNA data
#' @return A Seurat object which contains majority vote labels
#' @export
MapCellTypes <- function(sc_obj, reference, data_type = "scRNA") {
  message("Step 6: Reference-based cell type mapping...")
  anchors <- FindMappingAnchors(sc_obj, reference, data_type)
  sc_obj <- Seurat::MapQuery(anchorset = anchors,
                     query = sc_obj,
                     reference = reference,
                     refdata = "celltype.l2",
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     verbose = TRUE)
  sc_obj <- MajorityVote(sc_obj)
  return(sc_obj)
}
