% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_analysis.R
\name{run_downstream_analyses_RNA}
\alias{run_downstream_analyses_RNA}
\title{Perform downstream analyses (RNA)}
\usage{
run_downstream_analyses_RNA(
  sc_obj,
  reference_tissue,
  species,
  metadata_df,
  output_dir = getwd(),
  exit_with_code = FALSE,
  log_flag = log_flag
)
}
\arguments{
\item{sc_obj}{Seurat object containing cells for all samples}

\item{reference_tissue}{Reference tissue}

\item{species}{Species (human or mouse)}

\item{metadata_df}{Data frame containing metadata for samples. Rownames should be sample names and column names should be metadata attributes with two classes (e.g., condition: disease and control)}

\item{output_dir}{Path to directory where output will be saved. Defaults to working directory (\code{\link[=getwd]{getwd()}}). Directory will be created if it doesn't already exist.}

\item{exit_with_code}{Boolean flag to indicate whether we will terminate R session with exit code (via \code{\link[=quit]{quit()}}) if error occurs. If set to FALSE, we just use \code{\link[=stop]{stop()}}.}

\item{log_flag}{If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of \code{\link[=run_SPEEDI]{run_SPEEDI()}}.}
}
\value{
A list containing differential expression analyses
}
\description{
Perform downstream analyses (RNA)
}
\examples{
\dontrun{differential_expression_results <- RunDE_RNA(sc_obj = sc_obj, metadata_df = metadata_df)}
}
