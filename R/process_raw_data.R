#' Filter raw data for RNA
#'
#' @param all_sc_exp_matrices List of single cell expression matrices
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param record_doublets Boolean flag to indicate whether we will record doublets in the data (using the [scDblFinder] package). Possible choices are `TRUE` or `FALSE`.
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_file_path Path to log file (used to capture QC thresholds during parallel processing). Most likely only used in the context of [run_SPEEDI()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains filtered data from all input data
#' @examples
#' \dontrun{sc_obj <- FilterRawData_RNA(all_sc_exp_matrices)}
#' \dontrun{sc_obj <- FilterRawData_RNA(all_sc_exp_matrices,
#' species = "human", record_doublets = TRUE)}
#' @export
#' @importFrom foreach %dopar%
FilterRawData_RNA <- function(all_sc_exp_matrices, species = "human", record_doublets = FALSE, output_dir = getwd(), exit_with_code = FALSE, log_file_path = NULL, log_flag = FALSE) {
  exit_code <- -1
  sc_obj <- tryCatch(
    {
      # Normalize paths (in case user provides relative paths)
      output_dir <- normalize_dir_path(output_dir)
      species <- tolower(species)
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 3: Filtering out bad samples (RNA)", log_flag)
      print_SPEEDI(paste0("species is: ", species), log_flag)
      print_SPEEDI(paste0("record_doublets is: ", record_doublets), log_flag)
      print_SPEEDI("Creating Seurat object from matrices", log_flag)
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
      print_SPEEDI("Grabbing QC-related metrics", log_flag)
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
      # Print QC related output for user
      Create_QC_Output_Prefiltered_RNA(sc_obj, output_dir, log_flag)
      # Split up samples for individual processing
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

      print_SPEEDI("Beginning parallel processing of samples", log_flag)
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
            max_mt <- max(objects[[i]]$percent.mt)
          } else {
            max_mt <- kneedle::kneedle(graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$breaks[-1],
                                       graphics::hist(objects[[i]]$percent.mt, breaks=max(10, 0.5 * max(objects[[i]]$percent.mt)), plot=F)$counts)[1]
            max_mt <- max(max_mt, stats::quantile(objects[[i]]$percent.mt, .75))
          }
        } else { max_mt <- 0}
        # Subset current sample based on filtering criteria
        max_hb <- stats::quantile(objects[[i]]$percent.hb, .99)
        if (max_hb > 10) { max_hb <- 10 }
        if (max_hb < 2) { max_hb <- max(objects[[i]]$percent.hb) }

        max_rps <- stats::quantile(objects[[i]]$percent.rps, .99)
        if (max_rps < 2) { max_rps <- max(objects[[i]]$percent.rps) }

        max_rpl <- stats::quantile(objects[[i]]$percent.rpl, .99)
        if (max_rpl < 2) { max_rpl <- max(objects[[i]]$percent.rpl) }

        lower_bound_filtered_nFeature_cells <- length(objects[[i]]$nFeature_RNA[objects[[i]]$nFeature_RNA < lower_nF])
        upper_bound_filtered_nFeature_cells <- length(objects[[i]]$nFeature_RNA[objects[[i]]$nFeature_RNA >= stats::quantile(objects[[i]]$nFeature_RNA, .99)])
        upper_bound_filtered_percent.mt_cells <- length(objects[[i]]$percent.mt[objects[[i]]$percent.mt > max_mt])
        upper_bound_filtered_percent.rps_cells <- length(objects[[i]]$percent.rps[objects[[i]]$percent.rps > max_rps])
        upper_bound_filtered_percent.rpl_cells <- length(objects[[i]]$percent.rpl[objects[[i]]$percent.rpl > max_rpl])
        upper_bound_filtered_percent.hb_cells <- length(objects[[i]]$percent.hb[objects[[i]]$percent.hb > max_hb])
        object <- subset(x = objects[[i]],
                         subset = nFeature_RNA >= lower_nF &
                           nFeature_RNA < stats::quantile(objects[[i]]$nFeature_RNA, .99) &
                           percent.mt <= max_mt &
                           percent.rps <= max_rps &
                           percent.rpl <= max_rpl &
                           percent.hb <= max_hb)
        # Print info about QC thresholds for current sample to console
        # Note that this will not work in certain environments (e.g., RStudio) because of parallel processing
        message(paste0("QC Thresholds used for sample: ", current_sample_name))
        message(paste0("lower nFeature: ", lower_nF, " (", lower_bound_filtered_nFeature_cells, " cells failed QC)"))
        message(paste0("upper nFeature: ", stats::quantile(objects[[i]]$nFeature_RNA, .99), " (", upper_bound_filtered_nFeature_cells, " cells failed QC)"))
        message(paste0("max mt: ", max_mt, " (", upper_bound_filtered_percent.mt_cells, " cells failed QC)"))
        message(paste0("max rps: ", max_rps, " (", upper_bound_filtered_percent.rps_cells, " cells failed QC)"))
        message(paste0("max rpl: ", max_rpl, " (", upper_bound_filtered_percent.rpl_cells, " cells failed QC)"))
        message(paste0("max hb: ", max_hb, " (", upper_bound_filtered_percent.hb_cells, " cells failed QC)"))
        total_cells_filtered_out <- length(objects[[i]]$cell_name) - length(object$cell_name)
        message(paste0("Total cells filtered out: ", total_cells_filtered_out))
        message(paste0("Total cells remaining: ", length(object$cell_name)))
        # If we want to guarantee that our parameters are saved to a log file, we need a different strategy
        if(log_flag) {
          # Save sample parameters to a temporary text file
          current_sample_parameters <- paste0(current_sample_name, ",", lower_nF, ",", stats::quantile(objects[[i]]$nFeature_RNA, .99),
                                              ",", max_mt, ",", max_rps, ",", max_rpl,
                                              ",", max_hb, ",", lower_bound_filtered_nFeature_cells, ",", upper_bound_filtered_nFeature_cells,
                                              ",", upper_bound_filtered_percent.mt_cells, ",", upper_bound_filtered_percent.rps_cells, ",", upper_bound_filtered_percent.rpl_cells,
                                              ",", upper_bound_filtered_percent.hb_cells, ",", total_cells_filtered_out, ",", length(object$cell_name))
          sample_log_file_name <- paste0(log_file_path, "_", current_sample_name, ".QC.sample.txt")
          utils::write.table(current_sample_parameters, file = sample_log_file_name)
        }
        return(object)
      }
      if(log_flag) {
        # Read sample QC parameter files one at a time and write content into log
        sample_qc_file_paths <- list.files(dirname(log_file_path), pattern = paste0(basename(log_file_path), ".+.QC.sample.txt"), full.names = TRUE)
        for(sample_qc_file_path in sample_qc_file_paths) {
          sample_qc_stats <- strsplit(utils::read.table(sample_qc_file_path)$x, split = ",")[[1]]
          cat(paste0("\n", Sys.time(), ": QC Thresholds used for sample: ", sample_qc_stats[1]), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n\n", Sys.time(), ": lower nFeature: ", sample_qc_stats[2], " (", sample_qc_stats[8], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": upper nFeature: ", sample_qc_stats[3], " (", sample_qc_stats[9], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": max mt: ", sample_qc_stats[4], " (", sample_qc_stats[10], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": max rps: ", sample_qc_stats[5], " (", sample_qc_stats[11], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": max rpl: ", sample_qc_stats[6], " (", sample_qc_stats[12], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": max hb: ", sample_qc_stats[7], " (", sample_qc_stats[13], " cells failed QC)"), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": Total cells filtered out: ", sample_qc_stats[14]), file = paste0(log_file_path, ".log"), append = TRUE)
          cat(paste0("\n", Sys.time(), ": Total cells remaining: ", sample_qc_stats[15]), file = paste0(log_file_path, ".log"), append = TRUE)
          cat("\n", file = paste0(log_file_path, ".log"), append = TRUE)
          # Remove temporary file after we're done
          file.remove(sample_qc_file_path)
        }
        print_SPEEDI("\n", log_flag, silence_time = TRUE)
      }
      print_SPEEDI("Parallel processing complete", log_flag)

      print_SPEEDI(paste0("Filtered data has ", dim(sc_obj)[2], " barcodes and ", dim(sc_obj)[1], " transcripts."), log_flag)
      print_SPEEDI("Step 3: Complete", log_flag)
      return(sc_obj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running FilterRawData_RNA() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 16
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(sc_obj)
}

#' Filter raw data for ATAC
#'
#' @param proj ArchR project associated with data
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR project which contains filtered data
#' @examples
#' \dontrun{proj <- FilterRawData_ATAC(proj)}
#' @export
FilterRawData_ATAC <- function(proj, output_dir = getwd(), exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  proj <- tryCatch(
    {
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 3: Filtering out bad samples (ATAC)", log_flag)
      Create_QC_Output_Prefiltered_ATAC(proj, output_dir, log_flag)
      print_SPEEDI("Filtering out doublets and low quality cells (only keep cells which have TSS enrichment >= 12 and nucleosome ratio < 2)", log_flag)
      proj <- ArchR::filterDoublets(ArchRProj = proj)
      idxPass <- which(proj$TSSEnrichment >= 12 & proj$NucleosomeRatio < 2)
      cellsPass <- proj$cellNames[idxPass]
      proj <- proj[cellsPass, ]
      print_SPEEDI("Successfully filtered data", log_flag)
      print_SPEEDI("Step 3: Complete", log_flag)
      return(proj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running FilterRawData_ATAC() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 17
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(proj)
}

#' Process filtered data (RNA)
#'
#' @param sc_obj Seurat object containing cells for all samples
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]). Directory will be created if it doesn't already exist.
#' @param metadata_df Dataframe containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A Seurat object which contains processed data
#' @examples
#' \dontrun{sc_obj <- InitialProcessing_RNA(sc_obj)}
#' \dontrun{sc_obj <- InitialProcessing_RNA(sc_obj, species = "human")}
#' @export
InitialProcessing_RNA <- function(sc_obj, species = "human", output_dir = getwd(), metadata_df = NULL, exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  sc_obj <- tryCatch(
    {
      species <- tolower(species)
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 4: Processing raw data (RNA)", log_flag)
      print_SPEEDI(paste0("species is: ", species), log_flag)
      # Load cell cycle genes
      s.genes <- Seurat::cc.genes.updated.2019$s.genes
      g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
      m.s.genes <-  SPEEDI::cc.gene.updated.mouse$m.s.genes
      m.g2m.genes <-  SPEEDI::cc.gene.updated.mouse$m.g2m.genes
      # Perform initial SC transform (before cell cycle scoring)
      print_SPEEDI("Running SCTransform to normalize count data", log_flag)
      system.time(sc_obj <- Seurat::SCTransform(object = sc_obj,
                                                vst.flavor = "v2",
                                                vars.to.regress = c("percent.mt",
                                                                    "percent.rps",
                                                                    "percent.rpl",
                                                                    "percent.hb"),
                                                do.scale = TRUE,
                                                do.center = TRUE,
                                                return.only.var.genes = TRUE,
                                                seed.use = get_speedi_seed(),
                                                verbose = TRUE))
      # Perform cell cycle scoring
      print_SPEEDI("Loading cell cycle genes and performing cell cycle scoring", log_flag)
      if(species == "human") {
        sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
      } else {
        sc_obj <- Seurat::CellCycleScoring(object = sc_obj, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
      }
      sc_obj$CC.Difference <- sc_obj$S.Score - sc_obj$G2M.Score
      # Normalize count data using SCTransform
      print_SPEEDI("Running SCTransform to normalize count data (after performing cell cycle scoring)", log_flag)
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
                                                seed.use = get_speedi_seed(),
                                                verbose = TRUE))
      # Run PCA and UMAP to visualize data (prior to batch correction)
      print_SPEEDI("Running PCA and UMAP on normalized data", log_flag)
      sc_obj <- Seurat::RunPCA(sc_obj, npcs = 100, approx = T, seed.use = get_speedi_seed(), verbose = T)
      sc_obj <- Seurat::RunUMAP(sc_obj, reduction = "pca", dims = 1:100, seed.use = get_speedi_seed())
      sample_count <- length(unique(sc_obj$sample))
      cell_count <- length(sc_obj$cell_name)
      current_title <- paste0("RNA Data Before Integration\n(By Sample)\n(", sample_count, " Samples, ", cell_count, " Cells)")
      print_UMAP_RNA(sc_obj, file_name = "Before_Batch_Correction_RNA_UMAP_by_Sample.png",
                     group_by_category = "sample", output_dir = output_dir, title = current_title,
                     log_flag = log_flag)
      # Add metadata to samples
      if(!is.null(metadata_df)) {
        print_SPEEDI("Adding user metadata to samples", log_flag)
        for(i in 1:ncol(metadata_df)) {
          current_metadata_attribute <- colnames(metadata_df)[i]
          sample_metadata <- sc_obj$sample
          for(j in 1:nrow(metadata_df)) {
            sample_metadata <- gsub(rownames(metadata_df)[j], metadata_df[j,i], sample_metadata)
          }
          sc_obj <- Seurat::AddMetaData(object = sc_obj, metadata = sample_metadata, col.name = current_metadata_attribute)
        }
      }
      print_SPEEDI("Step 4: Complete", log_flag)
      return(sc_obj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running InitialProcessing_RNA() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 18
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(sc_obj)
}

#' Process filtered data (ATAC)
#'
#' @param proj ArchR project associated with data
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @param exit_with_code Boolean flag to indicate whether we will terminate R session with exit code (via [quit()]) if error occurs. If set to FALSE, we just use [stop()].
#' @return An ArchR project which contains processed data
#' @examples
#' \dontrun{proj <- InitialProcessing_ATAC(proj)}
#' @export
InitialProcessing_ATAC <- function(proj, output_dir = getwd(), exit_with_code = FALSE, log_flag = FALSE) {
  exit_code <- -1
  proj <- tryCatch(
    {
      # Normalize paths (in case user provides relative paths)
      output_dir <- normalize_dir_path(output_dir)
      print_SPEEDI("\n", log_flag, silence_time = TRUE)
      print_SPEEDI("Step 4: Processing raw data (ATAC)", log_flag)
      proj <- ArchR::addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",
                                     iterations = 2,
                                     force = TRUE,
                                     clusterParams = list(resolution = c(2), sampleCells = 10000, n.start = 30),
                                     varFeatures = 20000, dims = 1:30,
                                     saveIterations = TRUE)
      proj <- ArchR::addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
      proj <- ArchR::addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat",
                                 name = "Clusters", resolution = 2, knnAssign = 30,
                                 maxClusters = NULL, force = TRUE)
      print_SPEEDI("Printing UMAPs of filtered data", log_flag)
      num_cells <- length(proj$cellNames)
      num_samples <- length(unique(proj$Sample))
      sample_text <- ""
      if(num_samples == 1) {
        sample_text <- paste0("(1 Sample, ", num_cells, " Cells)")
      } else {
        sample_text <- paste0("(", num_samples, " Samples, ", num_cells, " Cells)")
      }
      p1 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("ATAC Data Before Integration\n(By Sample)\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), legend.text = ggplot2::element_text(size=10))
      ggplot2::ggsave(filename = paste0(output_dir, "Before_Batch_Correction_ATAC_UMAP_by_Sample.png"), plot = p1, device = "png", width = 8, height = 8, units = "in")
      p2 <- ArchR::plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "TSSEnrichment", embedding = "UMAP", force = TRUE, keepAxis = TRUE) +
        ggplot2::ggtitle(paste0("ATAC Data Before Integration\n(By TSS Enrichment)\n", sample_text)) + ggplot2::theme(plot.title = ggplot2::element_text(size=18), legend.key.size = ggplot2::unit(1, "cm"), legend.text = ggplot2::element_text(size=10))
      ggplot2::ggsave(filename = paste0(output_dir, "Before_Batch_Correction_ATAC_UMAP_by_TSSEnrichment.png"), plot = p2, device = "png", width = 8, height = 8, units = "in")
      # ArchR::plotPDF(p1,p2, name = "UMAPs_After_Initial_Processing_Plots", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
      print_SPEEDI("Step 4: Complete", log_flag)
      return(proj)
    },
    error = function(cond) {
      if(exit_code == -1) {
        print_SPEEDI("Error running InitialProcessing_ATAC() function", log_flag = log_flag)
        print_SPEEDI(cond, log_flag = log_flag)
        exit_code <- 19
      }
      quit_SPEEDI(exit_with_code = exit_with_code, exit_code = exit_code, log_flag = log_flag)
    }
  )
  gc()
  return(proj)
}
