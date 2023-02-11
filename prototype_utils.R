#------------------------------------------------
# Utility functions copied from Seurat Github
#  using OP data as an example
#
# Author: Yuan Wang
# Date:  08/31/2022
#------------------------------------------------

# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github('satijalab/seurat-wrappers')
# devtools::install_github("immunogenomics/lisi")
# devtools::install_github('theislab/kBET')
library(hdf5r)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratWrappers)
library(SingleCellExperiment)
library(plyr)
library(dplyr)
library(magrittr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(zoo)
library(ggfortify)
library(parallel)
library(purrr)
library(R.utils)
library(svMisc)
library(foreach)
library(doMC)
library(doParallel)
library(lisi)
library(outliers)
library(gmodels)
library(reticulate)
library(kBET)
library(irlba)
library(kneedle)
library(scCustomize)

# Check if the input matrix is dgCMatrix
#
# @param mat sparse matrix
# @return A dgCMatrix
#
RowSparseCheck <- function(mat) {
  if (!inherits(x = mat, what = "sparseMatrix")) {
    stop("Input should be sparse matrix")
  } else if (class(x = mat) != "dgCMatrix") {
    warning("Input matrix is converted to dgCMatrix.")
    mat <- as.sparse(x = mat)
  }
  return(mat)
}

# Calculate row variance of a sparse matrix
#
# @param mat sparse matrix
# @return A vector of row variance
#
RowVarSparse <- function(mat) {
  mat <- RowSparseCheck(mat = mat)
  output <- rowVars(mat)
  names(x = output) <- rownames(x = mat)
  return(output)
}


# Calculate row variance of a dense matrix
#
# @param mat dense matrix
# @return A vector of row variance
#
RowVar <- function(mat) {
  if (!inherits(x = mat, what = "matrix")) {
    stop("Input should be dense matrix")
  } else {
    output <- rowVars(mat)
    names(x = output) <- rownames(x = mat)
    return(output)
  }
}

# Scale a vector to range(0,1)
#
# @param x numeric vector
# @return A scaled vector ranging from 0 to 1
#
scale_zero_one <- function(x) {(x - min(x))/(max(x) - min(x))}


# Append a list to a list-of-lists
#
# @param lst, list
# @return A list of lists
#
lappend <- function (lst, ...){ c(lst, list(...))}

PrepDR <- function(object, features = NULL, slot = 'scale.data', verbose = TRUE) {
  if (length(x = VariableFeatures(object = object)) == 0 && is.null(x = features)) {
    stop("Variable features haven't been set. Run FindVariableFeatures() or provide a vector of feature names.")
  }
  data.use <- GetAssayData(object = object, slot = slot)
  if (nrow(x = data.use ) == 0 && slot == "scale.data") {
    stop("Data has not been scaled. Please run ScaleData and retry")
  }
  features <- features %||% VariableFeatures(object = object)
  features.keep <- unique(x = features[features %in% rownames(x = data.use)])
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have not been scaled (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep
  
  if (inherits(x = data.use, what = 'dgCMatrix')) {
    features.var <- RowVarSparse(mat = data.use[features, ])
  }
  else {
    features.var <- RowVar(mat = data.use[features, ])
  }
  features.keep <- features[features.var > 0]
  if (length(x = features.keep) < length(x = features)) {
    features.exclude <- setdiff(x = features, y = features.keep)
    if (verbose) {
      warning(paste0("The following ", length(x = features.exclude), " features requested have zero variance (running reduction without them): ", paste0(features.exclude, collapse = ", ")))
    }
  }
  features <- features.keep
  features <- features[!is.na(x = features)]
  data.use <- data.use[features, ]
  return(data.use)
}

PCA <- function(sc_obj) {
  # Fast PCA
  # sc_obj <- all.integrated.obj
  npcs <- min(30L, nrow(x = sc_obj))
  data.use <- PrepDR(object = sc_obj, features = NULL, verbose = TRUE)
  set.seed(SEED)
  system.time(pca_results <- prcomp_irlba(x = t(data.use), n = 30, retx = TRUE, center = TRUE, scale. = FALSE))
  feature.loadings <- pca_results$rotation
  sdev <- pca_results$sdev
  cell.embeddings <- pca_results$x / (pca_results$sdev[1:npcs] * sqrt(x = ncol(x = sc_obj) - 1))
  rownames(x = feature.loadings) <- rownames(x = sc_obj)
  colnames(x = feature.loadings) <- paste0("PC_", 1:npcs)
  rownames(x = cell.embeddings) <- colnames(x = sc_obj)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  
  total.variance <- sum(RowVar(mat = data.use))
  reduction.data <- CreateDimReducObject(
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = DefaultAssay(object = sc_obj),
    stdev = sdev,
    key = "PC_",
    misc = list(total.variance = total.variance))
  
  sc_obj[["prcomp"]] <- reduction.data
  
  return(sc_obj)
}

