LoadReference <- function(tissue, human) {
  if (human) {
    if (tissue == "Adipose") {
      SeuratData::InstallData("adiposeref")
      return(utils::data("adiposeref")) }

    if (tissue == "Bone Marrow") {
      SeuratData::InstallData("bonemarrowref")
      return(utils::data("bonemarrowref")) }

    if (tissue == "Fetus") {
      SeuratData::InstallData("fetusref")
      return(utils::data("fetusref")) }

    if (tissue == "Heart") {
      SeuratData::InstallData("heartref")
      return(utils::data("heartref")) }

    if (tissue == "Cortex") {
      SeuratData::InstallData("humancortexref")
      return(utils::data("humancortexref")) }

    if (tissue == "Kidney") {
      SeuratData::InstallData("kidneyref")
      return(utils::data("kidneyref")) }

    if (tissue == "Lung") {
      SeuratData::InstallData("lungref")
      return(utils::data("lungref")) }

    if (tissue == "Pancreas") {
      SeuratData::InstallData("pancreasref")
      return(utils::data("pancreasref")) }

    if (tissue == "PBMC") {
      reference_url <- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"
      reference <- SeuratDisk::LoadH5Seurat(RCurl::getURL(reference_url))
      return(reference) }

    if (tissue == "Tonsil") {
      SeuratData::InstallData("tonsilref")
      return(utils::data("tonsilref")) }
  }
  if (!human) {
    if (tissue == "Cortex") {
      SeuratData::InstallData("mousecortexref")
      return(utils::data("mousecortexref")) }
  }
}

FindMappingAnchors <- function(sc_obj, reference) {
  Seurat::DefaultAssay(sc_obj) <- "integrated"
  anchors <- Seurat::FindTransferAnchors(reference = reference,
                                 query = sc_obj,
                                 normalization.method = "SCT",
                                 recompute.residuals = T,
                                 reference.reduction = "spca")
  return(anchors)
}

MajorityVote <- function(sc_obj, current_resolution = 1.5) {
  associated_res_attribute <- paste0("integrated_snn_res.", current_resolution)
  message("Begin majority voting...")
  Seurat::DefaultAssay(sc_obj) <- "integrated"
  sc_obj <- Seurat::FindNeighbors(sc_obj, reduction = "pca", dims = 1:30)
  # TODO: Add code to find the best resolution (e.g., by using Clustree?)
  sc_obj <- Seurat::FindClusters(sc_obj, resolution = current_resolution)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)
  #integrated_snn_res_values <- sc_obj[[associated_res_attribute]]$

  cluster.dump <- as.numeric(levels(sc_obj$integrated_snn_res.1.5))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    print(i)
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(sc_obj$integrated_snn_res.1.5[cells]))
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
      cells <- names(sc_obj$integrated_snn_res.1.5[sc_obj$integrated_snn_res.1.5 == i])
      freq.table <- as.data.frame(table(sc_obj$predicted.id[cells]))
      levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(i)] <- as.vector(freq.table$Var1)[which.max(freq.table$Freq)]
    }
  }


  message("...End majority voting")
  return(sc_obj)
}

MapCellTypes <- function(sc_obj, reference) {
  message("Step 6: Reference-based cell type mapping...")
  anchors <- FindMappingAnchors(sc_obj, reference)
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
