% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{Read_RNA}
\alias{Read_RNA}
\title{Read in RNA data for processing}
\usage{
Read_RNA(
  input_dir = getwd(),
  sample_id_list = NULL,
  sample_file_paths = NULL,
  exit_with_code = FALSE,
  log_flag = FALSE
)
}
\arguments{
\item{input_dir}{Path to directory where input data are located. Defaults to working directory (\code{\link[=getwd]{getwd()}}).}

\item{sample_id_list}{Vector of sample names (optional - if not provided, will select all samples found recursively in \code{input_dir}).}

\item{sample_file_paths}{Vector of sample file paths (optional - if not provided, will select all samples found recursively in \code{input_dir}). If using Market Exchange (MEX) Format (matrix.mtx / barcodes.tsv / features.tsv or genes.tsv), please provide a full set of sample paths for only one type of file (e.g., \verb{"c("sample1/matrix.mtx", "sample2/matrix.mtx"}"). If this argument is used, \code{sample_id_list} is required and should be written in the same order as the sample file paths.}

\item{exit_with_code}{Boolean flag to indicate whether we will terminate R session with exit code (via \code{\link[=quit]{quit()}}) if error occurs. If set to FALSE, we just use \code{\link[=stop]{stop()}}.}

\item{log_flag}{If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of \code{\link[=run_SPEEDI]{run_SPEEDI()}}.}
}
\value{
A set of single cell expression matrices
}
\description{
Read in RNA data for processing
}
\examples{
\dontrun{all_sc_exp_matrices <- Read_RNA()}
\dontrun{all_sc_exp_matrices <- Read_RNA(input_dir = "~/input_data/",
sample_id_list = c("sample_1", "sample_2"))}
\dontrun{all_sc_exp_matrices <- Read_RNA(sample_id_list = c("sample_1", "sample_2"),
sample_file_paths = c("~/input_data/sample_1/filtered_feature_bc_matrix.h5",
"~/input_data/sample_2/filtered_feature_bc_matrix.h5"))}
}
