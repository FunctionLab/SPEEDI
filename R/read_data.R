#' Read in RNA data for processing
#'
#' @param data_path Path to directory where input data are located. Defaults to working directory ([getwd()]).
#' @param sample_id_list Vector of sample names (optional - if not provided, will select all samples found recursively in `data_path`).
#' @param sample_file_paths Vector of sample file paths (optional - if not provided, will select all samples found recursively in `data_path`). If using Market Exchange (MEX) Format (matrix.mtx / barcodes.tsv / features.tsv or genes.tsv), please provide a full set of sample paths for only one type of file (e.g., `"c("sample1/matrix.mtx", "sample2/matrix.mtx"`"). If this argument is used, `sample_id_list` is required and should be written in the same order as the sample file paths.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return A set of single cell expression matrices
#' @examples
#' \dontrun{all_sc_exp_matrices <- Read_RNA()}
#' \dontrun{all_sc_exp_matrices <- Read_RNA(data_path = "~/input_data/",
#' sample_id_list = c("sample_1", "sample_2"))}
#' \dontrun{all_sc_exp_matrices <- Read_RNA(sample_id_list = c("sample_1", "sample_2"),
#' sample_file_paths = c("~/input_data/sample_1/filtered_feature_bc_matrix.h5",
#' "~/input_data/sample_2/filtered_feature_bc_matrix.h5"))}
#' @export
#' @importFrom foreach %dopar%
Read_RNA <- function(data_path = getwd(), sample_id_list = NULL, sample_file_paths = NULL, log_flag = FALSE) {
  print_SPEEDI("Step 2 (RNA): Reading all samples", log_flag)
  # Normalize paths (in case user provides relative paths)
  data_path <- normalize_dir_path(data_path)
  # Print user parameters
  if(!is.null(data_path)) {
    print_SPEEDI(paste0("data_path is: ", data_path), log_flag)
  }
  if(!is.null(sample_id_list)) {
    print_SPEEDI("sample_id_list includes the following sample ids:", log_flag)
    print_SPEEDI(paste0(sample_id_list, collapse = ", "), log_flag)
  }
  if(!is.null(sample_file_paths)) {
    print_SPEEDI("sample_file_paths includes the following file paths:", log_flag)
    print_SPEEDI(paste0(sample_file_paths, collapse = ", "), log_flag)
  }
  if(is.null(sample_file_paths)) {
    # Make sure that data_path is fully expanded (aka replace ~ with full path to user's home dir)
    data_path <- path.expand(data_path)
    # First, remove "/" from end of data_path if it's provided (for use of list.files)
    last_char_of_data_path <- substr(data_path, nchar(data_path), nchar(data_path))
    if(last_char_of_data_path == "/") {
      data_path <- substr(data_path, 1, nchar(data_path) - 1)
    }
    # Second, look for all filtered_feature_bc_matrix .h5 files in data_path
    data_files <- list.files(path = data_path, pattern = "filtered_feature_bc_matrix\\.h5$", recursive = TRUE, full.names = TRUE)
    data_file_format <- "HDF5"
    if(length(data_files) == 0) {
      data_files <- list.files(path = data_path, pattern = "matrix.mtx", recursive = TRUE, full.names = TRUE)
      data_file_format <- "MEX"
    }
    # Finally, if the user did provide a sample_id_list, pick the subset of data files that have that sample ID in the path
    if(!is.null(sample_id_list)) {
      data_files <- data_files[grepl(paste(sample_id_list,collapse="|"), data_files)]
      positions <- sapply(sample_id_list, function(pattern) grep(pattern, data_files))
      sample_id_list <- sample_id_list[order(positions)]
    } else {
      if(length(data_files) > 1) {
        # Else, we are using all data files found above, but we still need to guess what the sample names are because of Cell Ranger's
        # structure for file output.
        # We use the following strategy: We look at the paths of all of our data files and split by "/"
        # Then, at each index, we see if the number of UNIQUE values is equivalent to the number of data files
        # If we find such an index, we assume that the sample names are stored there
        sample_id_list <- strsplit(data_files, "/")
        final_sample_id_list <- c()
        for(i in 1:length(sample_id_list[[1]])) {
          current_elements <- sapply(sample_id_list, "[[", i)
          if(length(unique(current_elements)) == length(data_files)) {
            final_sample_id_list <- current_elements
            break
          }
        }
        sample_id_list <- final_sample_id_list
      } else {
        # We can't use this strategy if users only provide one sample ID
        # Here, we'll just assume that the user's sample is the name of the folder right after data_path
        sample_id_list <- strsplit(data_files, paste0(data_path, "/"))
        sample_id_list <- sapply(sample_id_list , "[[", 2)
        sample_id_list <- strsplit(sample_id_list, "/")
        sample_id_list <- sapply(sample_id_list , "[[", 1)
      }
    }
  } else {
    data_files <- sample_file_paths
    positions <- sapply(sample_id_list, function(pattern) grep(pattern, data_files))
    sample_id_list <- sample_id_list[order(positions)]
  }
  print_SPEEDI(paste0("Total sample count is: ", length(sample_id_list)), log_flag)
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
  print_SPEEDI("Beginning parallel processing of samples", log_flag)
  # Dummy declaration to avoid check() complaining
  i <- 0
  all_sc_exp_matrices <- foreach::foreach(
    i = 1:length(data_files),
    .combine = 'cbind',
    .packages = c("Seurat", "base")
  ) %dopar% {
    # Read in data for current sample
    print(data_files[[i]])
    if(data_file_format == "HDF5") {
      sc_matrix <- Seurat::Read10X_h5(data_files[[i]])
    } else if(data_file_format == "MEX") {
      sc_matrix <- Seurat::Read10X(paste0(dirname(data_files[[i]]), "/"))
    }
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
  print_SPEEDI("Parallel processing complete", log_flag)
  print_SPEEDI(paste0("Raw data has ", dim(all_sc_exp_matrices)[2], " barcodes and ", dim(all_sc_exp_matrices)[1], " transcripts."), log_flag)
  print_SPEEDI("Step 2 (RNA): Complete", log_flag)
  gc()
  return(all_sc_exp_matrices)
}

#' Read in ATAC data for processing
#'
#' @param data_path Path to directory where input data are located. Defaults to working directory ([getwd()]).
#' @param output_dir Path to directory where output will be saved. Defaults to working directory ([getwd()]).
#' @param sample_id_list Vector of sample names (optional - if not provided, will select all samples found recursively in `data_path`).
#' @param sample_file_paths Vector of sample file paths (optional - if not provided, will select all samples found recursively in `data_path`). If used, `sample_id_list` is required.
#' @param species Species being analyzed. Possible choices are `"human"` or `"mouse"`.
#' @param log_flag If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of [run_SPEEDI()].
#' @return An ArchR project with associated Arrow files
#' @examples
#' \dontrun{proj <- Read_ATAC()}
#' \dontrun{proj <- Read_ATAC(data_path = "~/input_data/",
#' sample_id_list = c("sample_1", "sample_2"), species = "human")}
#' @export
#' @importFrom foreach %dopar%
Read_ATAC <- function(data_path = getwd(), output_dir = getwd(), sample_id_list = NULL, sample_file_paths = NULL, species = "human", log_flag = FALSE) {
  # Normalize paths (in case user provides relative paths)
  data_path <- normalize_dir_path(data_path)
  output_dir <- normalize_dir_path(output_dir)
  print_SPEEDI("Step 2 (ATAC): Reading all samples", log_flag)
  if(!is.null(data_path)) {
    print_SPEEDI(paste0("data_path is: ", data_path), log_flag)
  }
  if(!is.null(sample_id_list)) {
    print_SPEEDI(paste0("sample_id_list is: ", sample_id_list), log_flag)
  }
  if(!is.null(sample_file_paths)) {
    print_SPEEDI(paste0("sample_file_paths is: ", sample_file_paths), log_flag)
  }
  if(!is.null(sample_file_paths) & is.null(sample_id_list)) {
    print_SPEEDI("Error: You must provide a value for \"sample_id_list\" if you provide a value for \"sample_file_paths\".", log_flag)
    stop()
  }
  if(is.null(data_path) & is.null(sample_file_paths)) {
    print_SPEEDI("Error: You must provide a value for \"data_path\" if you do not provide a value for \"sample_file_paths\".", log_flag)
    stop()
  }
  if(is.null(sample_file_paths)) {
    # Make sure that data_path is fully expanded (aka replace ~ with full path to user's home dir)
    data_path <- path.expand(data_path)
    # First, remove "/" from end of data_path if it's provided (for use of list.files)
    last_char_of_data_path <- substr(data_path, nchar(data_path), nchar(data_path))
    if(last_char_of_data_path == "/") {
      data_path <- substr(data_path, 1, nchar(data_path) - 1)
    }
    # Second, look for all fragment.tsv.gz files in data_path
    data_files <- list.files(path = data_path, pattern = "fragments\\.tsv\\.gz$", recursive = TRUE, full.names = TRUE)
    # Finally, if the user did provide a sample_id_list, pick the subset of fragment files that have that sample ID in the path
    if(!is.null(sample_id_list)) {
      data_files <- data_files[grepl(paste(sample_id_list,collapse="|"), data_files)]
      positions <- sapply(sample_id_list, function(pattern) grep(pattern, data_files))
      sample_id_list <- sample_id_list[order(positions)]
    } else {
      # Else, we are using all data files found above, but we still need to guess what the sample names are because of Cell Ranger's
      # structure for file output.
      # We use the following strategy: We look at the paths of all of our data files and split by "/"
      # Then, at each index, we see if the number of UNIQUE values is equivalent to the number of data files
      # If we find such an index, we assume that the sample names are stored there
      sample_id_list <- strsplit(data_files, "/")
      final_sample_id_list <- c()
      for(i in 1:length(sample_id_list[[1]])) {
        current_elements <- sapply(sample_id_list, "[[", i)
        if(length(unique(current_elements)) == length(data_files)) {
          final_sample_id_list <- current_elements
          break
        }
      }
      sample_id_list <- final_sample_id_list
    }
  } else {
    data_files <- sample_file_paths
  }
  print_SPEEDI(paste0("Total sample count is: ", length(sample_id_list)), log_flag)
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
  # Set number of threads that ArchR will use
  ArchR::addArchRThreads(n.cores)
  # Set ArchR genome depending on species
  if(species == "human") {
    ArchR::addArchRGenome("hg38")
  } else {
    ArchR::addArchRGenome("mm10")
  }
  names(data_files) <- sample_id_list
  print_SPEEDI("Creating Arrow files (data files for ArchR)", log_flag)
  ArrowFiles <- ArchR::createArrowFiles(
    inputFiles = data_files,
    sampleNames = names(data_files),
    minTSS = 4,
    minFrags = 3000,
    maxFrags = 30000,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE
  )
  print_SPEEDI("Done creating Arrow files (data files for ArchR)", log_flag)
  # Calculate doublet scores
  print_SPEEDI("Calculating doublet scores", log_flag)
  doubScores <- ArchR::addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    LSIMethod = 1
  )
  print_SPEEDI("Done calculating doublet scores", log_flag)
  print_SPEEDI("Creating ArchR project", log_flag)
  proj <- ArchR::ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = TRUE
  )
  print_SPEEDI("Done creating ArchR project", log_flag)
  print_SPEEDI("Step 2 (ATAC): Complete", log_flag)
  gc()
  return(proj)
}
