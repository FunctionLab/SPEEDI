#------------------------------------------------
# Prototype API for TALOS
#
# Author: Yuan Wang
# Date:  08/08/2022
#------------------------------------------------

home_dir <- "/Genomics/argo/users/yuanwang/r04/yuanwang/DARPA/pipeline"
source(paste0(home_dir, "/prototype_utils.R"))

# Human CC Genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

# Mouse CC Genes
cc.gene.updated.mouse <- readRDS(paste0(home_dir, "/cc.gene.updated.mouse.rds"))
m.s.genes <- cc.gene.updated.mouse$m.s.genes
m.g2m.genes <- cc.gene.updated.mouse$m.g2m.genes


SEED <- 1824409L
set.seed(SEED)


Read_h5 <- function(data_path, sample_id_list) {
  message("Step 1: Reading all samples...")
  
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  if (n.cores > length(sample_id_list)) {
    n.cores <- length(sample_id_list)
  }
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallizing...")
  
  all_sc_exp_matrices <- foreach(
    i = 1:length(sample_id_list),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    library(hdf5r)
#     print(paste0(sample_id_list[[i]], "/filtered_feature_bc_matrix"))
#     sc_matrix <- Read10X(paste0(sample_id_list[[i]], "/filtered_feature_bc_matrix"))
    
    print(paste0(data_path, sample_id_list[[i]], "_filtered_feature_bc_matrix.h5"))
    sc_matrix <- Read10X_h5(paste0(data_path, sample_id_list[[i]], "_filtered_feature_bc_matrix.h5"))
    
    if (class(x = sc_matrix) == "list") {
      sc_exp_matrix <- sc_matrix$`Gene Expression`
    } else {
      sc_exp_matrix <- sc_matrix
    }
    if (grepl("_|\\.", i)) {
      prefix <- paste0(strsplit(sample_id_list[[i]], "_")[[1]][1], "#")
    } else {
      prefix <- paste0(sample_id_list[[i]], "#")
    }
    colnames(sc_exp_matrix) <- paste0(prefix, colnames(sc_exp_matrix))
    return(sc_exp_matrix)
  }
  
  message(paste0("Raw data has ", dim(all_sc_exp_matrices)[2], " barcodes and ", dim(all_sc_exp_matrices)[1], " transcripts."))
  return(all_sc_exp_matrices)
}

FilterRawData <- function(all_sc_exp_matrices, human) {
  message("Step 2: Filtering out bad samples...")
  
  sc_obj <- CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "unbias")
  
  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))
  
  if (human) {
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPS",
                                   col.name = "percent.rps") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPL",
                                   col.name = "percent.rpl") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^HB[A|B]",
                                   col.name = "percent.hb") 
    

  } else {
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^mt-",
                                   col.name = "percent.mt") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rps",
                                   col.name = "percent.rps") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rpl",
                                   col.name = "percent.rpl") 
    sc_obj <- PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Hb[a|b]",
                                   col.name = "percent.hb") 
  }
    
  objects <- SplitObject(sc_obj, split.by = "sample")

    
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  if (n.cores > length(objects)) {
    n.cores <- length(objects)
  }
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallizing...")
  
  sc_obj <- foreach(
    i = 1:length(objects),
    .combine = 'merge',
    .packages = c("Seurat", "base")
  ) %dopar% {  
 
#     lower_nF_dec <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
#                             hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts,
#                             decreasing = F)[1]
#     lower_nF_inc <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
#                             hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts,
#                             decreasing = T)[1]
#     lower_nF <- min(lower_nF_dec, lower_nF_inc)
      
    lower_nF <- kneedle(hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
                        hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts)[1]
    if (lower_nF > 1000) { lower_nF <- 1000 }
      
    if (max(objects[[i]]$percent.mt) > 0) {
        if (max(objects[[i]]$percent.mt) < 5) {
            max_mt <- quantile(objects[[i]]$percent.mt, .99)
        } else {
          max_mt <- kneedle(hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$breaks[-1],
                            hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$counts)[1]
          max_mt <- max(max_mt, quantile(objects[[i]]$percent.mt, .75))
        }
    } else { max_mt <- 0}
      
    max_hb <- quantile(objects[[i]]$percent.hb, .99)
    if (max_hb > 10) { max_hb <- 10 }
    
    object <- subset(x = objects[[i]],
                     subset = nFeature_RNA >= lower_nF &
                       nFeature_RNA < quantile(objects[[i]]$nFeature_RNA, .99) &
                       percent.mt <= max_mt &
                       percent.rps <= quantile(objects[[i]]$percent.rps, .99) &
                       percent.rpl <= quantile(objects[[i]]$percent.rpl, .99) &
                       percent.hb <= max_hb)

    return(object)
  }

  message(paste0("Filtered data has ", dim(sc_obj)[2], " barcodes and ", dim(sc_obj)[1], " transcripts."))
  return(sc_obj)
}

InitialProcessing <- function(sc_obj, human) {
  message("Step 3: Processing raw data...")
  if (human) {
      sc_obj <- CellCycleScoring(object = sc_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  } else {
       sc_obj <- CellCycleScoring(object = sc_obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  }
  sc_obj$CC.Difference <- sc_obj$S.Score - sc_obj$G2M.Score
  system.time(sc_obj <- SCTransform(object = sc_obj,
                                    vst.flavor = "v2",
                                    vars.to.regress = c("percent.mt",
                                                        "percent.rps",
                                                        "percent.rpl",
                                                        "percent.hb",
                                                        "CC.Difference"),
                                    do.scale = TRUE,
                                    do.center = TRUE,
                                    return.only.var.genes = TRUE,
                                    seed.use = SEED,
                                    verbose = TRUE))
#   sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 3000)
#   sc_obj <- ScaleData(sc_obj)
  sc_obj <- RunPCA(sc_obj, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:30, seed.use = SEED)
  return(sc_obj)
}

InferBatches <- function(sc_obj) {
  message("Step 4: Infer heterogeneous groups for integration...")
  
  sc_obj <- FindNeighbors(object = sc_obj, dims = 1:30)
  sc_obj <- FindClusters(object = sc_obj, resolution = 0.1, algorithm = 2, random.seed = SEED)
#   if (length(levels(sc_obj$seurat_clusters)) > 20) {
#       sc_obj <- FindClusters(object = sc_obj, resolution = 0.02, algorithm = 2, random.seed = SEED)
#   }
    
  
  # Use LISI metric to guess batch labels
  X <- sc_obj@reductions$umap@cell.embeddings
  meta_data <- data.frame(sc_obj$sample)
  colnames(meta_data) <- "batch"
  meta_data$cluster <- sc_obj$seurat_clusters
  lisi.res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(lisi.res) <- c("batch", "score", "cluster", "freq")
  clusters.interest <- names(table(sc_obj$seurat_clusters))[prop.table(table(sc_obj$seurat_clusters)) > 0.01]
  for (cluster in clusters.interest) { #levels(sc_obj$seurat_clusters)) {
    cells <- names(sc_obj$seurat_clusters[sc_obj$seurat_clusters == cluster])
    X.sub <- X[which(rownames(X) %in% cells),]
    meta_data.sub <- meta_data[which(rownames(meta_data) %in% cells),]
    res <- compute_lisi(X.sub, meta_data.sub, label_colnames = "batch")
    rownames(res) <- cells
    colnames(res) <- "score"
    res$batch <- meta_data.sub$batch
    agg.res <- aggregate(.~batch,data=res,mean)
    agg.res$cluster <- cluster
    agg.res$freq <- data.frame(table(res$batch))$Freq[which(data.frame(table(res$batch))$Var1 %in% agg.res$batch)]
    lisi.res <- rbind(lisi.res, agg.res)
  }
  
  p.values <- list()
  used.sample.dump <- c()
  batch.assign <- list()
  for ( i in clusters.interest) {
    #print(i)
    #lisi.res.sub <- lisi.res[lisi.res$score <= quantile(lisi.res$score, 1),]
    lisi.res.sub <- lisi.res[lisi.res$cluster == i,]
    if (max(lisi.res.sub$score) <= 1.1) {
      samples.of.batch <- lisi.res.sub$batch[1]
      if (!(samples.of.batch %in% used.sample.dump)) {
        batch.assign <- lappend(batch.assign, samples.of.batch)
      }
      used.sample.dump <- union(used.sample.dump, samples.of.batch)
    } else {
      lisi.res.sub$scaled.score <- scale_zero_one(lisi.res.sub$score * (lisi.res.sub$freq / sum(lisi.res.sub$freq)))
      lisi.res.sub <- lisi.res.sub[order(lisi.res.sub$scaled.score, decreasing = TRUE),]
      if (dim(lisi.res.sub)[1] > 30) {
        lisi.res.sub <- lisi.res.sub[1:30,]
      }
      lisi.res.sub$diff.scaled.score <- abs(c(diff(lisi.res.sub$scaled.score), 0))
        
      if (dim(lisi.res.sub)[1] >= 3) {
        p.values[[i]] <- dixon.test(lisi.res.sub$diff.scaled.score)$p.value[[1]]
      } else {
        p.values[[i]] <- 1
      }
    
      if (p.values[[i]] < 0.05) {
        max.index <- which.max(lisi.res.sub$diff.scaled.score)
        samples.of.batch <- lisi.res.sub$batch[1:max.index]
        
        if (any(samples.of.batch %in% used.sample.dump)) {
          if (!all(samples.of.batch %in% used.sample.dump)) {
            used.index <- which(samples.of.batch %in% used.sample.dump)
            samples.of.batch <- samples.of.batch[-used.index]
            if (length(samples.of.batch) > 0) {
                batch.assign <- lappend(batch.assign, samples.of.batch)
            }
          } else if (!list(samples.of.batch) %in% batch.assign) {
              if (length(samples.of.batch) == 1) {
                  batch.assign <- lappend(batch.assign, samples.of.batch)
              } else {
                  used.index <- which(samples.of.batch %in% unlist(batch.assign))
                  samples.of.batch <- samples.of.batch[-used.index]
                  if (length(samples.of.batch) > 0) {
                      batch.assign <- lappend(batch.assign, samples.of.batch)
                  }
              }
          }
        } else {
          batch.assign <- lappend(batch.assign, samples.of.batch)
        }
        used.sample.dump <- union(used.sample.dump, samples.of.batch)
      }
    }
  }
    
  batch <- as.factor(sc_obj$sample)
    
  if (length(batch.assign) > 0) {
      levels.batch <- levels(batch)
      for (i in 1:length(batch.assign)) { 
          levels.batch[which(levels(batch) %in% batch.assign[[i]])] <- i 
      }
      levels.batch[!levels.batch %in% c(1:length(batch.assign))] <- length(batch.assign)+1
      levels(batch) <- levels.batch
      sc_obj$batch <- as.character(batch)
  }
    
  else {
      message("No batch effect detected!")
      sc_obj$batch <- "No Batch"
  }

  print(unique(batch))
    
#   saveRDS(sc_obj, paste0(home_dir, "/unintegrated.object.RData"))
  return(sc_obj)
}

IntegrateByBatch <- function(sc_obj) {
  message("Step 5: Integrate samples based on inferred groups...")
  sc_obj_list <- SplitObject(sc_obj, split.by = "batch")
  
  
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  
  if (n.cores > length(sc_obj_list)) {
    n.cores <- length(sc_obj_list)
  }
  
  message(paste0("Number of cores: ", n.cores))
  
  registerDoMC(n.cores)
  message("Begin parallizing...")

  
  r <- foreach(
    i = 1:length(sc_obj_list),
    .combine = 'c',
    .packages = c("Seurat", "base")
  ) %dopar% {
    SEED <- 1824409L
      tmp <- SCTransform(object = sc_obj_list[[i]],
                       vst.flavor = "v2",
                       vars.to.regress = c("percent.mt",
                                           "percent.rps",
                                           "percent.rpl",
                                           "percent.hb",
                                           "CC.Difference"),
                       do.scale = TRUE,
                       do.center = TRUE,
                       return.only.var.genes = TRUE,
                       seed.use = SEED,
                       verbose = TRUE)
      tmp <- RunPCA(tmp, npcs = 30, approx = T, verbose = T, seed.use = SEED)
      return(tmp)
  }
  message(paste0(length(r), " samples transformed."))
  message("... Done parallizing")
  
  
  message("Select integration features...")
  features <- SelectIntegrationFeatures(object.list = r, nfeatures = 3000)
  r <- PrepSCTIntegration(object.list = r, anchor.features = features)
    
  message("Find integration anchors...")
  
  
#   anchors <- FindIntegrationAnchors(object.list = r,
#                                     normalization.method = "SCT",
#                                     anchor.features = features)

#   if (length(sc_obj_list) > 10) { 
#     message("...use reference-based integration...")
#     anchors <- FindIntegrationAnchors(object.list = r,
#                                       reference = 1,
#                                       normalization.method = "SCT",
#                                       anchor.features = features,
#                                       reduction = "rpca",
#                                       k.anchor = 10)
#   } else {
    anchors <- FindIntegrationAnchors(object.list = r,
                                      normalization.method = "SCT",
                                      anchor.features = features,
                                      reduction = "rpca",
                                      k.anchor = 10)
#    }
    
  message("Begin integration...")
#   integrated_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  integrated_obj <- IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT",
                                  k.weight = 100)
  DefaultAssay(integrated_obj) <- "integrated"
  
  rm(sc_obj_list)
  rm(features)
  rm(anchors)
  
  return(integrated_obj)
  # return(r)
}




VisualizeIntegration <- function(sc_obj) {
  set.seed(1824409L)
  sc_obj <- ScaleData(sc_obj, verbose = T)
#   sc_obj <- PCA(sc_obj)
  sc_obj <- RunPCA(sc_obj, npcs = 30, approx = T, verbose = T, seed.use = SEED)
  sc_obj <- RunUMAP(sc_obj, reduction = "pca", dims = 1:30, seed.use = SEED, return.model = T)
#   sc_obj <- RunUMAP(sc_obj, reduction = "prcomp", dims = 1:30, seed.use = SEED, return.model = T)
  DefaultAssay(sc_obj) <- "SCT"
  sc_obj <- PrepSCTFindMarkers(sc_obj)
  return(sc_obj)
}

LoadReference <- function(tissue, human) {
  if (human) {
    if (tissue == "Adipose") { 
        InstallData("adiposeref")
        return(data("adiposeref")) }
      
    if (tissue == "Bone Marrow") { 
        InstallData("bonemarrowref")
        return(data("bonemarrowref")) }
      
    if (tissue == "Fetus") { 
        InstallData("fetusref")
        return(data("fetusref")) }
      
    if (tissue == "Heart") { 
        InstallData("heartref")
        return(data("heartref")) }
      
    if (tissue == "Cortex") { 
        InstallData("humancortexref")
        return(data("humancortexref")) }
      
    if (tissue == "Kidney") { 
        InstallData("kidneyref")
        return(data("kidneyref")) }
      
    if (tissue == "Lung") { 
        InstallData("lungref")
        return(data("lungref")) }
      
    if (tissue == "Pancreas") { 
        InstallData("pancreasref")
        return(data("pancreasref")) }
      
    if (tissue == "PBMC") { 
        reference <- LoadH5Seurat(paste0(home_dir, "/reference/pbmc_multimodal.h5seurat"))
        return(reference) }
      
    if (tissue == "Tonsil") { 
        InstallData("tonsilref")
        return(data("tonsilref")) }
  }
  if (!human) {
    if (tissue == "Cortex") { 
        InstallData("mousecortexref")
        return(data("mousecortexref")) }
  }
}

FindMappingAnchors <- function(sc_obj, reference) {
  DefaultAssay(sc_obj) <- "integrated"
  anchors <- FindTransferAnchors(reference = reference,
                                 query = sc_obj,
                                 normalization.method = "SCT",
                                 recompute.residuals = T,
                                 reference.reduction = "spca")
  return(anchors)
}

MajorityVote <- function(sc_obj) {
  message("Begin majority voting...")
  DefaultAssay(sc_obj) <- "integrated"
  sc_obj <- FindNeighbors(sc_obj, reduction = "pca", dims = 1:30)
  sc_obj <- FindClusters(sc_obj, resolution = 1)
  sc_obj$predicted.id <- as.character(sc_obj$predicted.id)
    
  votes <- c()
  cluster.dump <- as.numeric(levels(sc_obj$integrated_snn_res.1))
  sc_obj$predicted_celltype_majority_vote <- sc_obj$seurat_clusters
  levels(sc_obj$predicted_celltype_majority_vote) <- as.character(levels(sc_obj$predicted_celltype_majority_vote))
  for (i in unique(sc_obj$predicted.id)) {
    cells <- names(sc_obj$predicted.id[sc_obj$predicted.id == i])
    freq.table <- as.data.frame(table(sc_obj$integrated_snn_res.1[cells]))
    freq.table <- freq.table[order(freq.table$Freq, decreasing = TRUE),]
    freq.table$diff <- abs(c(diff(freq.table$Freq), 0))
    p.values <- dixon.test(freq.table$diff)$p.value[[1]]
    max.index <- which.max(freq.table$diff)
    clusters <- as.numeric(as.character(freq.table$Var1[1:max.index]))
    levels(sc_obj$predicted_celltype_majority_vote)[levels(sc_obj$predicted_celltype_majority_vote) %in% as.character(clusters)] <- i
    cluster.dump <- cluster.dump[!cluster.dump %in% clusters]
  }
  
  if (length(cluster.dump) > 0) {
      for (i in cluster.dump) {
          cells <- names(sc_obj$integrated_snn_res.1[sc_obj$integrated_snn_res.1 == i])
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
  sc_obj <- MapQuery(anchorset = anchors,
                     query = sc_obj,
                     reference = reference,
                     refdata = "celltype.l2",
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     verbose = TRUE)
  sc_obj <- MajorityVote(sc_obj)
  return(sc_obj)
}




