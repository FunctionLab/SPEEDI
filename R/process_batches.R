#' Infer batches using LISI metric
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains labeled batches
#' @examples
#' \dontrun{sc_obj <- InferBatches(sc_obj)}
#' @export
#' @importFrom foreach %dopar%
InferBatches <- function(sc_obj, exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  sc_obj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 5: Inferring heterogeneous groups for integration", log_flag)
      # Find clusters in data (prior to batch correction)
      if ('lsi' %in% SeuratObject::Reductions(sc_obj)) {
        if (is.null(sc_obj@graphs$tileMatrix_snn)) {
          sc_obj <- Seurat::FindNeighbors(object = sc_obj, reduction = "lsi", dims = 1:100)
        } else {
          print_SPEEDI("Neighbors exist. Skipping construction of neighborhood graph...", log_flag)
        }
      } else {
        if (is.null(sc_obj@graphs$SCT_snn)) {
          sc_obj <- Seurat::FindNeighbors(object = sc_obj, reduction = "pca", dims = 1:100)
        } else {
          print_SPEEDI("Neighbors exist. Skipping construction of neighborhood graph...", log_flag)
        }
      }

      # Set up processing of different resolutions so it's parallel (max = 7 cores)
      if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
        n.cores <- as.numeric(parallel::detectCores())
      } else {
        n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
      }

      if (n.cores > 32) {
        n.cores <- 32
      }

      res <- seq(0.1, 1, by=0.01)

      print_SPEEDI(paste0("Number of cores: ", n.cores), log_flag)
      doParallel::registerDoParallel(n.cores)
      metrics_list <- foreach::foreach(
        i = res,
        .combine = 'c',
        .packages = c("Seurat", "bluster")
      ) %dopar% {
        tmp <- find_clusters_SPEEDI(sc_obj = sc_obj, resolution = i, method = "Louvain", log_flag = log_flag)
        dims <- 1:100
        clusters <- tmp$seurat_clusters
        if ('lsi' %in% SeuratObject::Reductions(sc_obj)) {
          sil.out <- bluster::approxSilhouette(Seurat::Embeddings(tmp@reductions$lsi)[, dims], clusters)
        } else {
          sil.out <- bluster::approxSilhouette(Seurat::Embeddings(tmp@reductions$pca)[, dims], clusters)
        }
        sil.score <- mean(sil.out$width)
        names(sil.score) <- i
        return(sil.score)
      }

      max.res <- as.numeric(as.character(names(which.max(metrics_list))))
      print_SPEEDI(paste0("Ideal Resolution = ", max.res), log_flag = log_flag)

      tmp <- find_clusters_SPEEDI(sc_obj = sc_obj, resolution = max.res, method = "Louvain", log_flag = log_flag)

      # Use LISI metric to guess batch labels
      X <- tmp@reductions$umap@cell.embeddings
      meta_data <- data.frame(tmp$sample)
      colnames(meta_data) <- "batch"
      meta_data$cluster <- tmp$seurat_clusters
      lisi.res <- data.frame(matrix(ncol = 4, nrow = 0))
      colnames(lisi.res) <- c("batch", "score", "cluster", "freq")
      clusters.interest <- names(table(tmp$seurat_clusters))[prop.table(table(tmp$seurat_clusters)) > 0.01]
      for (cluster in clusters.interest) {
        cells <- names(tmp$seurat_clusters[tmp$seurat_clusters == cluster])
        X.sub <- X[which(rownames(X) %in% cells),]
        meta_data.sub <- meta_data[which(rownames(meta_data) %in% cells),]
        res <- lisi::compute_lisi(X.sub, meta_data.sub, label_colnames = "batch")
        rownames(res) <- cells
        colnames(res) <- "score"

        res$score <- res$score + apply(X.sub, 1, dist)
        res$score <- res$score + max(abs(res$score)) + 0.01

        res$batch <- meta_data.sub$batch
        agg.res <- stats::aggregate(.~batch,data=res,mean)
        agg.res$cluster <- cluster
        agg.res$freq <- data.frame(table(res$batch))$Freq[which(data.frame(table(res$batch))$Var1 %in% agg.res$batch)]
        lisi.res <- rbind(lisi.res, agg.res)
      }

      lisi.res <- lisi.res[lisi.res$freq > 10,]

      p.values <- list()
      used.sample.dump <- c()
      batch.assign <- list()
      for ( i in clusters.interest) {
        lisi.res.sub <- lisi.res[lisi.res$cluster == i,]
        lisi.res.sub$batch_count <- as.numeric(as.character(table(tmp$sample)[lisi.res.sub$batch]))

        lisi.res.sub$scaled.score <- mapply(
          geometric.mean,
          (max(lisi.res.sub$score) / lisi.res.sub$score),
          (lisi.res.sub$freq / lisi.res.sub$batch_count))

        if (max(lisi.res.sub$score) <= 1.01) {
          lisi.res.sub <- lisi.res.sub[order(lisi.res.sub$score, decreasing = TRUE),]
          samples.of.batch <- lisi.res.sub$batch[1]

          if (!(list(samples.of.batch) %in% batch.assign)) {
            batch.assign <- lappend(batch.assign, samples.of.batch)
          }
          used.sample.dump <- union(used.sample.dump, samples.of.batch)
        } else {
          lisi.res.sub <- lisi.res.sub[order(lisi.res.sub$scaled.score, decreasing = TRUE),]

          if (dim(lisi.res.sub)[1] > 30) {
            lisi.res.sub <- lisi.res.sub[1:30,]
          }

          lisi.res.sub$diff.scaled.score <- abs(c(diff(lisi.res.sub$scaled.score), 0))

          if (dim(lisi.res.sub)[1] >= 3 & sum(lisi.res.sub$diff.scaled.score) != 0) {
            p.values[[i]] <- outliers::dixon.test(lisi.res.sub$diff.scaled.score)$p.value[[1]]
          } else {
            p.values[[i]] <- 1
          }

          if (p.values[[i]] < 0.05 & dim(lisi.res.sub)[1] >= 3) {
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

      batch <- as.factor(tmp$sample)

      if (length(batch.assign) > 0) {
        levels.batch <- levels(batch)
        for (i in 1:length(batch.assign)) {
          levels.batch[which(levels(batch) %in% batch.assign[[i]])] <- i
        }
        levels.batch[!levels.batch %in% c(1:length(batch.assign))] <- length(batch.assign)+1
        levels(batch) <- levels.batch
        tmp$batch <- as.character(batch)
      } else {
        print_SPEEDI("No batch effect detected!", log_flag)
        tmp$batch <- "No Batch"
      }
      print_SPEEDI(paste0("Total batches detected: ", length(unique(tmp$batch))), log_flag)
      print_SPEEDI("Step 5: Complete", log_flag)
      return(tmp)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running InferBatches() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 20
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(sc_obj)
}

#' Integrate batches (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains integrated data
#' @examples
#' \dontrun{sc_obj <- IntegrateByBatch_RNA(sc_obj)}
#' @export
#' @importFrom foreach %dopar%
IntegrateByBatch_RNA <- function(sc_obj, exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  integrated_obj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 6: Integrating samples based on inferred groups (RNA)", log_flag)
      # If we only have one batch, we don't need to integrate by batch, so we exit the function
      if(length(unique(sc_obj$batch)) == 1) {
        print_SPEEDI("Only one batch was found, so we don't need to integrate batches. Exiting IntegrateByBatch!", log_flag)
        return(sc_obj)
      }
      DefaultAssay(sc_obj) <- "RNA"
      sc_obj <- JoinLayers(sc_obj)
      sc_obj[["RNA"]] <- split(sc_obj[["RNA"]], f = sc_obj$batch)
      regress_vars <- c("percent.mt", "percent.rps", "percent.rpl", "percent.hb", "CC.Difference")
      regress_vars <- regress_vars[which(regress_vars %in% colnames(sc_obj_list[[i]][[]]))]
      if (length(regress_vars) == 0) { regress_vars <- NULL }
      sc_obj <- SCTransform(object = sc_obj,
                                      vst.flavor = "v2",
                                      vars.to.regress = regress_vars,
                                      do.scale = TRUE,
                                      do.center = TRUE,
                                      return.only.var.genes = TRUE,
                                      seed.use = get_speedi_seed(),
                                      verbose = TRUE)
      sc_obj <- RunPCA(sc_obj, npcs = 100, approx = T, verbose = T, seed.use = get_speedi_seed())
      sc_obj <- IntegrateLayers(
        object = sc_obj, method = RPCAIntegration,
        orig.reduction = "pca", new.reduction = "integrated.rpca",
        normalization.method = "SCT",
        dims.to.integrate = 100,
        k.weight = 90,
        verbose = TRUE
      )

      print_SPEEDI("Finished integration", log_flag)
      print_SPEEDI("Step 6: Complete", log_flag)
      return(integrated_obj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running IntegrateByBatch_RNA() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 21
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(integrated_obj)
}

#' Integrate batches (ATAC)
#'
#' @param proj ArchR project containing cells for all samples
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR object which contains integrated data
#' @examples
#' \dontrun{proj <- IntegrateByBatch_ATAC(proj)}
#' @export
IntegrateByBatch_ATAC <- function(proj, output_dir = getwd(), exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  proj <- tryCatch(
    {
      # Normalize paths (in case user provides relative paths)
      output_dir <- normalize_dir_path(output_dir)
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Preparing ATAC samples for batch inference", log_flag)
      tile_sce <- ArchR::getMatrixFromProject(proj, useMatrix='TileMatrix', binarize = TRUE)
      tile_reduc <- ArchR::getReducedDims(ArchRProj = proj, reducedDims = "IterativeLSI", returnMatrix = TRUE)
      tile_reduc <- tile_reduc[match(colnames(tile_sce), rownames(tile_reduc)),]
      for (i in colnames(SummarizedExperiment::colData(tile_sce))) {
        SummarizedExperiment::colData(tile_sce)[[i]] <- as.vector(SummarizedExperiment::colData(tile_sce)[[i]])
      }
      rownames(tile_sce) <- paste0(as.character(SummarizedExperiment::rowData(tile_sce)$seqnames),
                                   '-',
                                   as.character(SummarizedExperiment::rowData(tile_sce)$start))
      tile_seurat <- Seurat::CreateSeuratObject(SummarizedExperiment::assays(tile_sce)$TileMatrix[, rownames(tile_reduc)],
                                                project = "peaks",
                                                assay = "tileMatrix")
      # Doesn't currently work, but I don't think it's necessary
      # tile_seurat <- Seurat::AddMetaData(tile_seurat, data.frame(t(SummarizedExperiment::colData(tile_sce))))
      # LSI
      cell.embeddings <- tile_reduc
      feature.loadings <- matrix()
      assay <- "tileMatrix"
      sdev <- 0
      reduction.key <- "LSI_"
      reduction.data <- Seurat::CreateDimReducObject(
        embeddings = cell.embeddings,
        loadings = feature.loadings,
        assay = assay,
        stdev = sdev,
        key = reduction.key,
        misc = list()
      )
      tile_seurat@reductions$lsi <- reduction.data
      # UMAP
      tile_umap <- ArchR::getEmbedding(ArchRProj = proj, embedding = "UMAP", returnDF = TRUE)
      cell.embeddings_UMAP <- as.matrix(tile_umap)
      cell.embeddings_UMAP <- cell.embeddings_UMAP[match(rownames(cell.embeddings), rownames(cell.embeddings_UMAP)),]
      feature.loadings <- matrix()
      assay <- "tileMatrix"
      sdev <- 0
      reduction.key <- "UMAP_"
      reduction.data <- Seurat::CreateDimReducObject(
        embeddings = cell.embeddings_UMAP,
        loadings = feature.loadings,
        assay = assay,
        stdev = sdev,
        key = reduction.key,
        misc = list()
      )
      tile_seurat@reductions$umap <- reduction.data
      tile_seurat$sample <- proj$Sample
      # Step 4: Inferring batches
      tile_seurat <- InferBatches(tile_seurat, exit_with_code = exit_with_code, log_flag = log_flag)
      proj$Batch <- tile_seurat$batch
      # If we only have one batch, we don't need to integrate by batch, so we exit the function
      if(length(unique(proj$Batch)) == 1) {
        print_SPEEDI("Only one batch was found, so we don't need to integrate batches. Skipping Step 6 (integration of samples) and Step 7 Part 1 (final processing of integrated data).", log_flag)
      } else {
        print_SPEEDI("\n", log_flag, silence_time = TRUE)
        print_SPEEDI("Step 6: Integrating samples based on inferred groups", log_flag)
        print_SPEEDI("Beginning integration", log_flag)
        proj <- ArchR::addHarmony(ArchRProj = proj, reducedDims = "IterativeLSI",
                                  dimsToUse = 2:30, scaleDims = TRUE,
                                  corCutOff = 0.75, groupBy = "Batch", force = TRUE)
        print_SPEEDI("Step 6: Complete", log_flag)
        print_SPEEDI("\n", log_flag, silence_time = TRUE)
        print_SPEEDI("Step 7 (Part 1): Final processing of integrated data (ATAC)", log_flag)
        proj <- ArchR::addUMAP(ArchRProj = proj, reducedDims = "Harmony", force = TRUE)
        proj <- ArchR::addClusters(input = proj, reducedDims = "Harmony", method = "Seurat",
                                   name = "Harmony_clusters", resolution = 2, knnAssign = 30,
                                   maxClusters = NULL, force = TRUE)
        print_SPEEDI("Step 7: Complete", log_flag)
      }
      print_SPEEDI("Step 7 (Part 2): Clustering with Seurat and printing final UMAPs (ATAC)", log_flag)
      if(length(unique(proj$Batch)) == 1) {
        reducedDims_param <- "IterativeLSI"
      } else {
        reducedDims_param <- "Harmony"
      }
      tile_reduc <- ArchR::getReducedDims(ArchRProj = proj, reducedDims = reducedDims_param, returnMatrix = TRUE)
      tmp <- matrix(stats::rnorm(nrow(tile_reduc) * 3, 10), ncol = nrow(tile_reduc), nrow = 3)
      colnames(tmp) <- rownames(tile_reduc)
      rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
      obj <- Seurat::CreateSeuratObject(tmp, project='scATAC', min.cells=0, min.features=0)
      obj[[reducedDims_param]] <- Seurat::CreateDimReducObject(embeddings=tile_reduc, key=paste0(reducedDims_param, "_"), assay='RNA')
      obj <- Seurat::FindNeighbors(obj, reduction = reducedDims_param, dims = 1:29)
      obj <- find_clusters_SPEEDI(obj, resolution = 2, method = "Leiden", log_flag = log_flag)
      obj <- Seurat::RunUMAP(obj, reduction = reducedDims_param, dims = 1:29, seed.use = get_speedi_seed())
      proj <- ArchR::addCellColData(
        ArchRProj = proj,
        cells = names(obj$seurat_clusters),
        data = as.character(obj$seurat_clusters),
        name = "seurat_clusters",
        force = TRUE)
      rm(obj)
      # Grab info for titles of plots
      num_cells <- length(proj$cellNames)
      num_samples <- length(unique(proj$Sample))
      sample_text <- ""
      if(num_samples == 1) {
        sample_text <- paste0("(1 Sample, ", num_cells, " Cells)")
      } else {
        sample_text <- paste0("(", num_samples, " Samples, ", num_cells, " Cells)")
      }
      # Plot integrated UMAPs
      p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "seurat_clusters", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Clusters)\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), legend.text = ggplot2::element_text(size=10))
      ggplot2::ggsave(filename = paste0(output_dir, "Final_ATAC_UMAP_by_Clusters.png"), plot = p1, device = "png", width = 8, height = 8, units = "in")
      p2 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("ATAC Data Integration\n(By Sample)\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), legend.text = ggplot2::element_text(size=10))
      ggplot2::ggsave(filename = paste0(output_dir, "Final_ATAC_UMAP_by_Sample.png"), plot = p2, device = "png", width = 8, height = 8, units = "in")
      p3 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("ATAC Data Integration\n(By TSS Enrichment)\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), legend.key.size = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size=10))
      ggplot2::ggsave(filename = paste0(output_dir, "Final_ATAC_UMAP_by_TSSEnrichment.png"), plot = p3, device = "png", width = 8, height = 8, units = "in")
      # ArchR::plotPDF(p1,p2,p3, name = "UMAP_Final_Integrated_Plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
      return(proj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running IntegrateByBatch_ATAC() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 22
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(proj)
}

#' Visualize integration and prepare SCT for finding markers
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param resolution Resolution used for clustering
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object with SCT markers and visualizations
#' @examples
#' \dontrun{sc_obj <- VisualizeIntegration(sc_obj)}
#' @export
VisualizeIntegration <- function(sc_obj, output_dir = getwd(), resolution = 2, exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  sc_obj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 7: Final processing of integrated data (RNA)", log_flag)
      print_SPEEDI("Preparing integrated data for FindMarkers", log_flag)
      sc_obj <- Seurat::PrepSCTFindMarkers(sc_obj)
      print_SPEEDI("Finding clusters and printing UMAPs of integrated data", log_flag)
      if(length(unique(sc_obj$batch)) != 1) {
        reduction <- "integrated.rpca"
        umap_reduction_name <- "umap.rpca"
      } else {
        reduction <- "pca"
        umap_reduction_name <- "umap"
      }
      sc_obj <- Seurat::FindNeighbors(object = sc_obj, reduction = reduction, dims = 1:100)
      sc_obj <- find_clusters_SPEEDI(sc_obj = sc_obj, resolution = resolution, method = "Louvain", log_flag = log_flag)
      sc_obj <- Seurat::RunUMAP(sc_obj, reduction = reduction, dims = 1:100, seed.use = get_speedi_seed(), return.model = T, reduction.name = umap_reduction_name)
      sample_count <- length(unique(sc_obj$sample))
      cell_count <- length(sc_obj$cell_name)
      # Plot by cluster
      current_title <- paste0("RNA Data Integration\n(By Cluster)\n(", sample_count, " Samples, ", cell_count, " Cells)")
      print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Cluster.png",
                     group_by_category = "seurat_clusters", output_dir = output_dir, title = current_title,
                     log_flag = log_flag)
      # Plot by sample
      current_title <- paste0("RNA Data Integration\n(By Sample)\n(", sample_count, " Samples, ", cell_count, " Cells)")
      print_UMAP_RNA(sc_obj, file_name = "Final_RNA_UMAP_by_Sample.png",
                     group_by_category = "sample", output_dir = output_dir, title = current_title,
                     log_flag = log_flag)
      print_SPEEDI("Step 7: Complete", log_flag)
      return(sc_obj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running VisualizeIntegration() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 23
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(sc_obj)
}
