# We use the packages below, but not with direct calls in our functions.
# Thus, to avoid check() complaining about us not using these packages in
# its import check, we add some direct calls in the (never used) function below
ignore_unused_imports <- function() {
  hdf5r::as_hex
  Matrix::Arith
  glmGamPoi::glm_gp
  RCurl::AUTH_ANY
}
