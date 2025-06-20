# Environment to store global package variables
the <- new.env(parent = emptyenv())

# URL for PBMC reference - if the link becomes defunct and we haven't updated it yet,
# update the variable below (either by the setter below or by directly editing this file)
the$pbmc_reference_url <- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"

# URL for checking which genes are on HumanBase (and getting associated ENTREZ IDs)
the$gene_check_url <- 'https://hb.flatironinstitute.org/api/genes/search/multi?'

# URL for submitting functional module discovery job to HumanBase
the$fmd_submission_url <- 'https://hb.flatironinstitute.org/api/integrations/community/?integration='

# List of possible references
the$references <- c("adipose", "bone_marrow", "fetus",
                               "heart", "cortex", "kidney",
                               "lung", "pancreas", "pbmc", "pbmc_full", "tonsil", "custom", "none")

# List of Azimuth references
the$azimuth_references <- c("adiposeref", "bonemarrowref", "fetusref",
                    "heartref", "humancortexref", "kidneyref",
                    "lungref", "pancreasref", "pbmcref", "tonsilref", "mousecortexref")

# SPEEDI random seed
the$speedi_seed <- 1253


#' Get PBMC reference URL
#' @export
get_pbmc_reference_url <- function() {
  the$pbmc_reference_url
}

#' Set PBMC reference URL
#' @param new_url the new PBMC reference url
#' @export
set_pbmc_reference_url <- function(new_url = "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat") {
  old <- the$pbmc_reference_url
  the$pbmc_reference_url <- new_url
  invisible(old)
}

#' Get HumanBase gene check URL
#' @export
get_gene_check_url <- function() {
  the$gene_check_url
}

#' Get FMD submission URL
#' @export
get_fmd_submission_url <- function() {
  the$fmd_submission_url
}

#' Get possible reference options
#' @export
get_references <- function() {
  the$references
}

#' Get possible Azimuth reference options
#' @export
get_azimuth_references <- function() {
  the$azimuth_references
}

#' Get SPEEDI random seed
#' @export
get_speedi_seed <- function() {
  the$speedi_seed
}

#' Set SPEEDI random seed
#' @param new_seed the new SPEEDI random seed
#' @export
set_speedi_seed <- function(new_seed = 1253) {
  old <- the$speedi_seed
  the$speedi_seed <- new_seed
  invisible(old)
}

# We use the packages below, but not with direct calls in our functions.
# Thus, to avoid check() complaining about us not using these packages in
# its import check, we add some direct calls in the (never used) function below
ignore_unused_imports <- function() {
  BiocManager::version
  Biostrings::AAString
  chromVAR::addGCBias
  ComplexHeatmap::add_heatmap
  DESeq2::counts
  dplyr::across
  edgeR::addPriorCount
  GenomicRanges::absoluteRanges
  glmGamPoi::glm_gp
  hdf5r::as_hex
  hexbin::BTC
  limma::alias2Symbol
  Matrix::Arith
  motifmatchr::match_motifs
  RCurl::AUTH_ANY
  rhdf5::H5Aclose
  Rsamtools::applyPileups
  SeuratObject::Assays
}
