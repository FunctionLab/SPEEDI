#' Read in data
#'
#' @param data_path path to where data is located
#' @param sample_id_list list of sample names
#' @return A set of single cell expression matrices
#' @export
#' @importFrom foreach %dopar%
Read_h5 <- function(data_path, sample_id_list) {
  message("Step 1: Reading all samples...")
  # Make reading data parallel
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(parallel::detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }

  if (n.cores > length(sample_id_list)) {
    n.cores <- length(sample_id_list)
  }
  message(paste0("Number of cores: ", n.cores))
  doParallel::registerDoParallel(n.cores)
  message("Begin parallelizing")
  i <- 0
  all_sc_exp_matrices <- foreach::foreach(
    i = 1:length(sample_id_list),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    print(paste0(data_path, sample_id_list[[i]], "/outs/filtered_feature_bc_matrix.h5"))
    sc_matrix <- Seurat::Read10X_h5(paste0(data_path, sample_id_list[[i]], "/outs/filtered_feature_bc_matrix.h5"))
    if (inherits(x = sc_matrix, what = 'list')) {
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
