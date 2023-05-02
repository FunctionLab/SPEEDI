#' Read in h5 data for processing
#'
#' @param data_path path to where data is located (optional - if not provided, will assume the data path is the current working directory)
#' @param sample_id_list list of sample names (optional - if not provided, will select all samples found recursively in data_path)
#' @param log_flag if set to TRUE, we previously set up a log file where certain output will be written (e.g., parameters)
#' @return A set of single cell expression matrices
#' @export
#' @importFrom foreach %dopar%
Read_h5 <- function(data_path = getwd(), sample_id_list = NULL, log_flag = FALSE) {
  print_SPEEDI("Step 1: Reading all samples", log_flag)
  print_SPEEDI(paste0("data_path is: ", data_path), log_flag)
  if(!is.null(sample_id_list)) {
    print_SPEEDI(paste0("sample_id_list is: ", sample_id_list), log_flag)
  }
  # Make sure that data_path is fully expanded (aka replace ~ with full path to user's home dir)
  data_path <- path.expand(data_path)
  # First, remove "/" from end of data_path if it's provided (for use of list.files)
  last_char_of_data_path <- substr(data_path, nchar(data_path), nchar(data_path))
  if(last_char_of_data_path == "/") {
    data_path <- substr(data_path, 1, nchar(data_path) - 1)
  }
  # Second, look for all filtered_feature_bc_matrix .h5 files in data_path
  data_files <- list.files(path = data_path, pattern = "filtered_feature_bc_matrix\\.h5$", recursive = TRUE, full.names = TRUE)
  # Finally, if the user did provide a sample_id_list, pick the subset of .h5 files that have that sample ID in the path
  if(!is.null(sample_id_list)) {
    data_files <- data_files[grepl(paste(sample_id_list,collapse="|"), data_files)]
  } else {
    # Else, we are using all data files found above, but we still need to guess what the sample names are because of Cell Ranger's
    # structure for file output.
    # Our current approach assumes that sample names are the directories right after data_path.
    # Is there a better way of doing this?
    sample_id_list <- strsplit(data_files, paste0(data_path, "/"))
    sample_id_list <- sapply(sample_id_list , "[[", 2)
    sample_id_list <- strsplit(sample_id_list, "/")
    sample_id_list <- sapply(sample_id_list , "[[", 1)
  }
  # Set up reading of data so it's parallel (max cores == number of samples)
  if (Sys.getenv("SLURM_NTASKS_PER_NODE") == "") {
    n.cores <- as.numeric(parallel::detectCores())
  } else {
    n.cores <- as.numeric(Sys.getenv("SLURM_NTASKS_PER_NODE"))
  }
  if (n.cores > length(data_files)) {
    n.cores <- length(data_files)
  }
  print_SPEEDI(paste0("Number of cores: ", n.cores))
  doParallel::registerDoParallel(n.cores)
  print_SPEEDI("Begin parallelizing")
  # Dummy declaration to avoid check() complaining
  i <- 0
  all_sc_exp_matrices <- foreach::foreach(
    i = 1:length(data_files),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    # Read in data for current sample
    print(data_files[[i]])
    sc_matrix <- Seurat::Read10X_h5(data_files[[i]])
    # If our resulting data structure is a list, then we have some multiome data (we just want the gene expression)
    if (inherits(x = sc_matrix, what = 'list')) {
      sc_exp_matrix <- sc_matrix$`Gene Expression`
    } else {
      sc_exp_matrix <- sc_matrix
    }
    # Add sample names as prefix to cell names
    if (grepl("_|\\.", i)) {
      prefix <- paste0(strsplit(sample_id_list[[i]], "_")[[1]][1], "#")
    } else {
      prefix <- paste0(sample_id_list[[i]], "#")
    }
    colnames(sc_exp_matrix) <- paste0(prefix, colnames(sc_exp_matrix))
    return(sc_exp_matrix)
  }
  print_SPEEDI("\n", log_flag, silence_time = TRUE)
  print_SPEEDI(paste0("Raw data has ", dim(all_sc_exp_matrices)[2], " barcodes and ", dim(all_sc_exp_matrices)[1], " transcripts."), log_flag)
  print_SPEEDI("Step 1: Complete", log_flag)
  gc()
  return(all_sc_exp_matrices)
}
