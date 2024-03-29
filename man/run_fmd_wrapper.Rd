% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downstream_analysis.R
\name{run_fmd_wrapper}
\alias{run_fmd_wrapper}
\title{Wrapper for running FMD job in the context of the SPEEDI wrapper function}
\usage{
run_fmd_wrapper(
  gene_list,
  network,
  RNA_output_dir,
  cell_type,
  metadata_attribute,
  fc_flag,
  log_flag = FALSE
)
}
\arguments{
\item{gene_list}{Input list of genes}

\item{network}{Background network}

\item{RNA_output_dir}{Output directory (RNA)}

\item{cell_type}{Cell type}

\item{metadata_attribute}{Metadata attribute}

\item{fc_flag}{Flag to indicate direction of fold change ("high" or "low")}

\item{log_flag}{If set to TRUE, record certain output (e.g., parameters) to a previously set up log file. Most likely only used in the context of \code{\link[=run_SPEEDI]{run_SPEEDI()}}.}
}
\value{
A list containing FMD results
}
\description{
Wrapper for running FMD job in the context of the SPEEDI wrapper function
}
\examples{
\dontrun{hb_networks <- grab_hb_networks(current_cell_type = "B naive",
reference_tissue = "pbmc")}
}
