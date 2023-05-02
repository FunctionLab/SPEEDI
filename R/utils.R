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

#' Print to console as well as log file (if it's present)
#' @param current_message Message to print
#' @param log_flag boolean to indicate whether we're also printing to log file
#' @return TRUE
print_SPEEDI <- function(current_message, log_flag = FALSE) {
  message(current_message)
  if(log_flag) {
    logr::log_print(current_message, console = FALSE, hide_notes = TRUE, blank_after = FALSE)
  }
  return(TRUE)
}
