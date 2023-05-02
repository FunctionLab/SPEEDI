# Environment to store global package variables
the <- new.env(parent = emptyenv())

# URL for PBMC reference - if the link becomes defunct and we haven't updated it yet,
# update the variable below (either by the setter below or by directly editing this file)
the$pbmc_reference_url <- "https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat"

# List of possible SeuratData references
the$seuratdata_references <- c("adiposeref", "bonemarrowref", "fetusref",
                               "heartref", "humancortexref", "kidneyref",
                               "lungref", "pancreasref", "pbmcref", "tonsilref", "mousecortexref")

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

#' Get possible SeuratData references
#' @export
get_seuratdata_references <- function() {
  the$seuratdata_references
}

# We use the packages below, but not with direct calls in our functions.
# Thus, to avoid check() complaining about us not using these packages in
# its import check, we add some direct calls in the (never used) function below
ignore_unused_imports <- function() {
  hdf5r::as_hex
  Matrix::Arith
  glmGamPoi::glm_gp
  RCurl::AUTH_ANY
}
