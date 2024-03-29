% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_raw_data.R
\name{InitialProcessing_ATAC}
\alias{InitialProcessing_ATAC}
\title{Process filtered data (ATAC)}
\usage{
InitialProcessing_ATAC(
  proj,
  output_dir = getwd(),
  exit_with_code = FALSE,
  log_flag = FALSE
)
}
\arguments{
\item{proj}{ArchR project associated with data}

\item{output_dir}{Path to directory where output will be saved. Defaults to working directory (\code{\link[=getwd]{getwd()}}).}

\item{exit_with_code}{Boolean flag to indicate whether we will terminate R session with exit code (via \code{\link[=quit]{quit()}}) if error occurs. If set to FALSE, we just use \code{\link[=stop]{stop()}}.}

\item{log_flag}{If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of \code{\link[=run_SPEEDI]{run_SPEEDI()}}.}
}
\value{
An ArchR project which contains processed data
}
\description{
Process filtered data (ATAC)
}
\examples{
\dontrun{proj <- InitialProcessing_ATAC(proj)}
}
