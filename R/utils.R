#' Scale a vector to range(0,1)
#'
#' @param x Numeric vector
#' @return A scaled vector ranging from 0 to 1
scale_zero_one <- function(x) {(x - min(x))/(max(x) - min(x))}

#' Append a list to a list-of-lists
#'
#' @param lst List
#' @param ... Additional lists
#' @return A list of lists
lappend <- function (lst, ...){ c(lst, list(...))}
