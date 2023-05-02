#' Infer batches using LISI metric
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A Seurat object which contains labeled batches
#' @export
InferBatches <- function(sc_obj, log_flag = FALSE) {
  print_SPEEDI("\n", log_flag)
  print_SPEEDI("Step 4: Inferring heterogeneous groups for integration", log_flag)
  print_SPEEDI(paste0("log_flag is: ", log_flag), log_flag)
  # Find clusters in data (prior to batch correction)
  sc_obj <- Seurat::FindNeighbors(object = sc_obj, dims = 1:30)
  sc_obj <- Seurat::FindClusters(object = sc_obj, resolution = 0.1, algorithm = 2)
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
    res <- lisi::compute_lisi(X.sub, meta_data.sub, label_colnames = "batch")
    rownames(res) <- cells
    colnames(res) <- "score"
    res$batch <- meta_data.sub$batch
    agg.res <- stats::aggregate(.~batch,data=res,mean)
    agg.res$cluster <- cluster
    agg.res$freq <- data.frame(table(res$batch))$Freq[which(data.frame(table(res$batch))$Var1 %in% agg.res$batch)]
    lisi.res <- rbind(lisi.res, agg.res)
  }

  p.values <- list()
  used.sample.dump <- c()
  batch.assign <- list()
  for ( i in clusters.interest) {
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
        p.values[[i]] <- outliers::dixon.test(lisi.res.sub$diff.scaled.score)$p.value[[1]]
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
  } else {
    print_SPEEDI("No batch effect detected!", log_flag)
    sc_obj$batch <- "No Batch"
  }

  print_SPEEDI(paste0("Total batches detected: ", length(unique(batch))), log_flag)
  return(sc_obj)
}

#' Integrate batches
#'
#' @param sc_obj Seurat object
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A Seurat object which contains integrated data
#' @export
#' @importFrom foreach %dopar%
IntegrateByBatch <- function(sc_obj, log_flag = FALSE) {
  print_SPEEDI("\n", log_flag)
  print_SPEEDI("Step 5: Integrating samples based on inferred groups", log_flag)
  print_SPEEDI(paste0("log_flag is: ", log_flag), log_flag)
  sc_obj_list <- Seurat::SplitObject(sc_obj, split.by = "batch")
  # If we only have one batch, we don't need to integrate by batch, so we exit the function
  if(length(sc_obj_list) == 1) {
    print_SPEEDI("Only one batch was found, so we don't need to integrate batches. Exiting IntegrateByBatch!", log_flag)
    return(sc_obj)
  }
  # Set up reading of data so it's parallel (max cores == number of samples)
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(parallel::detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(sc_obj_list)) {
    n.cores <- length(sc_obj_list)
  }

  print_SPEEDI(paste0("Number of cores: ", n.cores), log_flag)

  doParallel::registerDoParallel(n.cores)
  print_SPEEDI("Begin parallelizing...", log_flag)
  # Dummy declaration to avoid check() complaining
  i <- 0
  r <- foreach::foreach(
    i = 1:length(sc_obj_list),
    .combine = 'c',
    .packages = c("Seurat", "base")
  ) %dopar% {
    # Normalize counts within batch
    tmp <- Seurat::SCTransform(object = sc_obj_list[[i]],
                       vst.flavor = "v2",
                       vars.to.regress = c("percent.mt",
                                           "percent.rps",
                                           "percent.rpl",
                                           "percent.hb",
                                           "CC.Difference"),
                       do.scale = TRUE,
                       do.center = TRUE,
                       return.only.var.genes = TRUE,
                       verbose = TRUE)
    tmp <- Seurat::RunPCA(tmp, npcs = 30, approx = T, verbose = T)
    return(tmp)
  }
  print_SPEEDI(paste0(length(r), " samples transformed."), log_flag)
  print_SPEEDI("... Done parallelizing...", log_flag)


  print_SPEEDI("Selecting integration features", log_flag)
  features <- Seurat::SelectIntegrationFeatures(object.list = r, nfeatures = 3000)
  r <- Seurat::PrepSCTIntegration(object.list = r, anchor.features = features)

  print_SPEEDI("Finding integration anchors", log_flag)

  anchors <- Seurat::FindIntegrationAnchors(object.list = r,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    reduction = "rpca",
                                    k.anchor = 10)

  print_SPEEDI("Beginning integration", log_flag)
  integrated_obj <- Seurat::IntegrateData(anchorset = anchors,
                                  normalization.method = "SCT",
                                  k.weight = 100)
  Seurat::DefaultAssay(integrated_obj) <- "integrated"

  rm(sc_obj_list)
  rm(features)
  rm(anchors)

  return(integrated_obj)
}


#' Visualize integration and prepare SCT for finding markers
#'
#' @param sc_obj Seurat object
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A Seurat object with SCT markers and visualizations
#' @export
VisualizeIntegration <- function(sc_obj, log_flag = FALSE) {
  print_SPEEDI("\n", log_flag)
  print_SPEEDI("Step 6: Scaling integrated data, creating UMAP of integration and prepping data for FindMarkers", log_flag)
  print_SPEEDI(paste0("log_flag is: ", log_flag), log_flag)
  sc_obj <- Seurat::ScaleData(sc_obj, verbose = T)
  sc_obj <- Seurat::RunPCA(sc_obj, npcs = 30, approx = T, verbose = T)
  sc_obj <- Seurat::RunUMAP(sc_obj, reduction = "pca", dims = 1:30, return.model = T)
  Seurat::DefaultAssay(sc_obj) <- "SCT"
  sc_obj <- Seurat::PrepSCTFindMarkers(sc_obj)
  print_SPEEDI("Step 6 completed", log_flag)
  return(sc_obj)
}

