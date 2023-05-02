#' Filter raw data
#'
#' @param all_sc_exp_matrices list of single cell expression matrices
#' @param species flag to indicate whether we're processing human or mouse data
#' @param record_doublets flag to indicate whether we're recording doublet info in the data (using scDblFinder)
#' @param log_file_path path to log file (used to capture QC thresholds during parallel processing)
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A Seurat object which contains filtered data from all input data
#' @export
#' @importFrom foreach %dopar%
FilterRawData <- function(all_sc_exp_matrices, species = "human", record_doublets = FALSE, log_file_path = NULL, log_flag = FALSE) {
  species <- tolower(species)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 2: Filtering out bad samples", log_flag)
  print_SPEEDI(paste0("species is: ", species), log_flag)
  print_SPEEDI(paste0("record_doublets is: ", record_doublets), log_flag)
  sc_obj <- Seurat::CreateSeuratObject(counts = all_sc_exp_matrices,
                               assay = "RNA",
                               min.cells = 3,
                               min.features = 3,
                               project = "unbias")
  # Label each cell with its sample name in the sample metadata column
  # Also, label each cell with its cell name in the cell_names metadata column
  sc_obj$sample <- as.vector(sapply(strsplit(colnames(sc_obj), "#"), "[", 1))
  cell_names <- rownames(sc_obj@meta.data)
  sc_obj <- Seurat::AddMetaData(sc_obj, metadata = cell_names, col.name = "cell_name")
  # Dummy declaration to avoid check() complaining
  scDblFinder.class <- NULL
  # If requested, record doublet info in samples
  if(record_doublets) {
    print_SPEEDI("Recording doublet info", log_flag)
    sc_obj <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(sc_obj), samples = "sample"))
    # See distribution of doublets in each sample
    doublet_sc_obj <- subset(x = sc_obj, subset = scDblFinder.class %in% "doublet")
    print_SPEEDI("Number of doublets in each sample:", log_flag)
    print_SPEEDI(table(doublet_sc_obj$sample), log_flag)
    rm(doublet_sc_obj)
  }
  # Grab QC-related metrics (regular expression is different for human vs. mouse)
  if (species == "human") {
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

  # Set up reading of data so it's parallel (max cores == number of samples)
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(parallel::detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(objects)) {
    n.cores <- length(objects)
  }

  print_SPEEDI(paste0("Number of cores: ", n.cores), log_flag)
  doParallel::registerDoParallel(n.cores)
  print_SPEEDI("Begin parallelizing", log_flag)
  # Dummy declarations to avoid check() complaining
  i <- 0
  nFeature_RNA <- percent.mt <- percent.rps <- percent.rpl <- percent.hb <- NULL
  sample_parameters <- list()
  sc_obj <- foreach::foreach(
    i = 1:length(objects),
    .combine = 'merge',
    .packages = c("Seurat", "base")
  ) %dopar% {
    current_sample_name <- unique(objects[[i]]$sample)
    # Automatically detect lower bound for nFeature using kneedle
    lower_nF <- kneedle::kneedle(graphics::hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$breaks[-1],
                                 graphics::hist(objects[[i]]$nFeature_RNA, breaks=100, plot=F)$counts)[1]
    if (lower_nF > 1000) { lower_nF <- 1000 }
    # Automatically detect upper bound for percent.mt using kneedle
    if (max(objects[[i]]$percent.mt) > 0) {
      if (max(objects[[i]]$percent.mt) < 5) {
        max_mt <- stats::quantile(objects[[i]]$percent.mt, .99)
      } else {
        max_mt <- kneedle::kneedle(graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$breaks[-1],
                                   graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$counts)[1]
        max_mt <- max(max_mt, stats::quantile(objects[[i]]$percent.mt, .75))
      }
    } else { max_mt <- 0}
    # Subset current sample based on filtering criteria
    max_hb <- stats::quantile(objects[[i]]$percent.hb, .99)
    if (max_hb > 10) { max_hb <- 10 }

    object <- subset(x = objects[[i]],
                      subset = nFeature_RNA >= lower_nF &
                       nFeature_RNA < stats::quantile(objects[[i]]$nFeature_RNA, .99) &
                       percent.mt <= max_mt &
                       percent.rps <= stats::quantile(objects[[i]]$percent.rps, .99) &
                       percent.rpl <= stats::quantile(objects[[i]]$percent.rpl, .99) &
                       percent.hb <= max_hb)
    # Print info about QC thresholds for current sample to console
    # Note that this will not work in certain environments (e.g., RStudio) because of parallel processing
    message(paste0("QC Thresholds used for sample: ", current_sample_name))
    message(paste0("lower nFeature: ", lower_nF), log_flag)
    message(paste0("upper nFeature: ", stats::quantile(objects[[i]]$nFeature_RNA, .99)))
    message(paste0("max mt: ", max_mt), log_flag)
    message(paste0("max rps: ", stats::quantile(objects[[i]]$percent.rps, .99)))
    message(paste0("max rpl: ", stats::quantile(objects[[i]]$percent.rpl, .99)))
    message(paste0("max hb: ", max_hb))
    # If we want to guarantee that our parameters are saved to a log file, we need a different strategy
    if(log_flag) {
      # Save sample parameters to a temporary text file
      current_sample_parameters <- paste0(current_sample_name, ",", lower_nF, ",", stats::quantile(objects[[i]]$nFeature_RNA, .99),
                                          ",", max_mt, ",", stats::quantile(objects[[i]]$percent.rps, .99), ",", stats::quantile(objects[[i]]$percent.rpl, .99),
                                          ",", max_hb)
      sample_log_file_name <- paste0(log_file_path, "_", current_sample_name, ".QC.sample.txt")
      write.table(current_sample_parameters, file = sample_log_file_name)
    }
    return(object)
  }
  if(log_flag) {
    # Read sample QC parameter files one at a time and write content into log
    sample_qc_file_paths <- list.files(dirname(log_file_path), pattern = paste0(basename(log_file_path), ".+.QC.sample.txt"), full.names = TRUE)
    for(sample_qc_file_path in sample_qc_file_paths) {
      sample_qc_stats <- strsplit(read.table(sample_qc_file_path)$x, split = ",")[[1]]
      cat(paste0("\n", Sys.time(), ": QC Thresholds used for sample: ", sample_qc_stats[1]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n\n", Sys.time(), ": lower nFeature: ", sample_qc_stats[2]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n", Sys.time(), ": upper nFeature: ", sample_qc_stats[3]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n", Sys.time(), ": max mt: ", sample_qc_stats[4]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n", Sys.time(), ": max rps: ", sample_qc_stats[5]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n", Sys.time(), ": max rpl: ", sample_qc_stats[6]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat(paste0("\n", Sys.time(), ": max hb: ", sample_qc_stats[7]), file = paste0(log_file_path, ".log"), append = TRUE)
      cat("\n", file = paste0(log_file_path, ".log"), append = TRUE)
      # Remove temporary file after we're done
      file.remove(sample_qc_file_path)
    }
    print_SPEEDI("\n", log_flag, silence_time = TRUE)
  }
  print_SPEEDI(paste0("Filtered data has ", dim(sc_obj)[2], " barcodes and ", dim(sc_obj)[1], " transcripts."), log_flag)
  print_SPEEDI("Step 2: Complete", log_flag)
  gc()
  return(sc_obj)
}

#' Process filtered data
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param species flag to indicate whether we're processing human or mouse data
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A Seurat object which contains processed data
#' @export
InitialProcessing <- function(sc_obj, species = "human", log_flag = FALSE) {
  species <- tolower(species)
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI("Step 3: Processing raw data", log_flag)
  print_SPEEDI(paste0("species is: ", species), log_flag)
  # Load cell cycle genes and perform cell cycle scoring
  s.genes <- Seurat::cc.genes.updated.2019$s.genes
  g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
  m.s.genes <-  SPEEDI::cc.gene.updated.mouse$m.s.genes
  m.g2m.genes <-  SPEEDI::cc.gene.updated.mouse$m.g2m.genes
  if (species == "human") {
    sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  } else {
    sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  }
  sc_obj$CC.Difference <- sc_obj$S.Score - sc_obj$G2M.Score
  # Normalize count data using SCTransform
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
  # Run PCA and UMAP to visualize data (prior to batch correction)
  # TODO: Print plot?
  sc_obj <- Seurat::RunPCA(sc_obj, npcs = 30, approx = T, verbose = T)
  sc_obj <- Seurat::RunUMAP(sc_obj, reduction = "pca", dims = 1:30)
  print_SPEEDI("Step 3: Complete", log_flag)
  gc()
  return(sc_obj)
}
