# We use the imports below, but not with direct calls in our functions.
# Thus, to avoid check() complaining about us not using imports, we add
# some (never used) direct calls for the imports below.
ignore_unused_imports <- function() {
  hdf5r::as_hex
  Matrix::Arith
  glmGamPoi::glm_gp
  RCurl::AUTH_ANY
}
