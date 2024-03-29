% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_data.R
\name{Read_ATAC}
\alias{Read_ATAC}
\title{Read in ATAC data for processing}
\usage{
Read_ATAC(
  input_dir = getwd(),
  output_dir = getwd(),
  sample_id_list = NULL,
  sample_file_paths = NULL,
  species = "human",
  exit_with_code = FALSE,
  log_flag = FALSE
)
}
\arguments{
\item{input_dir}{Path to directory where input data are located. Defaults to working directory (\code{\link[=getwd]{getwd()}}).}

\item{output_dir}{Path to directory where output will be saved. Defaults to working directory (\code{\link[=getwd]{getwd()}}).}

\item{sample_id_list}{Vector of sample names (optional - if not provided, will select all samples found recursively in \code{input_dir}).}

\item{sample_file_paths}{Vector of sample file paths (optional - if not provided, will select all samples found recursively in \code{input_dir}). If used, \code{sample_id_list} is required.}

\item{species}{Species being analyzed. Possible choices are \code{"human"} or \code{"mouse"}.}

\item{exit_with_code}{Boolean flag to indicate whether we will terminate R session with exit code (via \code{\link[=quit]{quit()}}) if error occurs. If set to FALSE, we just use \code{\link[=stop]{stop()}}.}

\item{log_flag}{If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of \code{\link[=run_SPEEDI]{run_SPEEDI()}}.}
}
\value{
An ArchR project with associated Arrow files
}
\description{
Read in ATAC data for processing
}
\examples{
\dontrun{proj <- Read_ATAC()}
\dontrun{proj <- Read_ATAC(input_dir = "~/input_data/",
sample_id_list = c("sample_1", "sample_2"), species = "human")}
}
