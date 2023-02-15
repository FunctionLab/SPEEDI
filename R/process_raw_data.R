#' Filter raw data
#'
#' @param all_sc_exp_matrices list of single cell expression matrices
#' @param human flag to indicate whether we're processing human or mouse data
#' @return A Seurat object which contains filtered data from all input data
#' @export
#' @importFrom foreach %dopar%
FilterRawData <- function(all_sc_exp_matrices, human) {
  message("Step 2: Filtering out bad samples...")
  sc_obj <- Seurat::CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "unbias")

  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))

  if (human) {
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPS",
                                   col.name = "percent.rps")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RPL",
                                   col.name = "percent.rpl")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^HB[A|B]",
                                   col.name = "percent.hb")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^RP[SL]",
                                   col.name = "percent.rp")


  } else {
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^mt-",
                                   col.name = "percent.mt")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rps",
                                   col.name = "percent.rps")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rpl",
                                   col.name = "percent.rpl")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Hb[a|b]",
                                   col.name = "percent.hb")
    sc_obj <- Seurat::PercentageFeatureSet(object = sc_obj,
                                   pattern = "^Rp[sl]",
                                   col.name = "percent.rp")
  }

  objects <- Seurat::SplitObject(sc_obj, split.by = "sample")


  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(parallel::detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(objects)) {
    n.cores <- length(objects)
  }

  message(paste0("Number of cores: ", n.cores))
  doParallel::registerDoParallel(n.cores)
  message("Begin parallelizing...")
  i <- 0
  nFeature_RNA <- percent.mt <- percent.rps <- percent.rpl <- percent.hb <- NULL
  sc_obj <- foreach::foreach(
    i = 1:length(objects),
    .combine = 'merge',
    .packages = c("Seurat", "base")
  ) %dopar% {
    current_sample_name <- unique(objects[[i]]$sample)
    lower_nF <- kneedle::kneedle(graphics::hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
                                 graphics::hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts)[1]
    if (lower_nF > 1000) { lower_nF <- 1000 }

    if (max(objects[[i]]$percent.mt) > 0) {
      if (max(objects[[i]]$percent.mt) < 5) {
        max_mt <- stats::quantile(objects[[i]]$percent.mt, .99)
      } else {
        max_mt <- kneedle::kneedle(graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$breaks[-1],
                                   graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$counts)[1]
        max_mt <- max(max_mt, stats::quantile(objects[[i]]$percent.mt, .75))
      }
    } else { max_mt <- 0}

    max_hb <- stats::quantile(objects[[i]]$percent.hb, .99)
    if (max_hb > 10) { max_hb <- 10 }

    object <- subset(x = objects[[i]],
                      subset = nFeature_RNA >= lower_nF &
                       nFeature_RNA < stats::quantile(objects[[i]]$nFeature_RNA, .99) &
                       percent.mt <= max_mt &
                       percent.rps <= stats::quantile(objects[[i]]$percent.rps, .99) &
                       percent.rpl <= stats::quantile(objects[[i]]$percent.rpl, .99) &
                       percent.hb <= max_hb)
    return(object)
  }

  message(paste0("Filtered data has ", dim(sc_obj)[2], " barcodes and ", dim(sc_obj)[1], " transcripts."))
  return(sc_obj)
}

#' Process filtered data
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param human flag to indicate whether we're processing human or mouse data
#' @return A Seurat object which contains processed data
#' @export
InitialProcessing <- function(sc_obj, human) {
  message("Step 3: Processing raw data...")
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  # TODO: Add this file to data dir?
  #cc.gene.updated.mouse <- readRDS(paste0(home_dir, "/cc.gene.updated.mouse.rds"))
  #m.s.genes <-  Seurat::cc.gene.updated.mouse$m.s.genes
  #m.g2m.genes <-  Seurat::cc.gene.updated.mouse$m.g2m.genes
  if (human) {
    sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  } else {
    sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  }
  sc_obj$CC.Difference <- sc_obj$S.Score - sc_obj$G2M.Score
  system.time(sc_obj <- Seurat::SCTransform(object = sc_obj,
                                    vst.flavor = "v2",
                                    vars.to.regress = c("percent.mt",
                                                        "percent.rps",
                                                        "percent.rpl",
                                                        "percent.hb",
                                                        "CC.Difference"),
                                    do.scale = TRUE,
                                    do.center = TRUE,
                                    return.only.var.genes = TRUE,
                                    verbose = TRUE))
  sc_obj <- Seurat::RunPCA(sc_obj, npcs = 30, approx = T, verbose = T)
  sc_obj <- Seurat::RunUMAP(sc_obj, reduction = "pca", dims = 1:30)
  return(sc_obj)
}
