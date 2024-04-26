####################################################################
# There is a bug in the original logr log_warning, so we make a fix here.
####################################################################

# Finally got this working
#' @title Logs a warning
#' @description Writes a warning message to the log. Warning will be written
#' both to the log at the point the function is called, and also written to the
#' message file.  This function is used internally.
#' @param msg The message to log.
#' @seealso \code{\link{log_error}} to write an error message to the log.
#' @export
#' @examples
#' library(logr)
#'
#' # Create temp file location
#' tmp <- file.path(tempdir(), "test.log")
#'
#' # Open log
#' lf <- log_open(tmp)
#'
#' # Send warning message to log
#' log_warning("Here is a warning")
#'
#' # Close log
#' log_close()
#'
#' # View results
#' writeLines(readLines(lf))
#'
log_warning <- function(msg = NULL) {

  logr:::update_status()


  if (e$log_status == "open") {


    # print("Warning Handler")
    if (!is.null(msg)) {
      msg1 <- paste("Warning:", msg)
      log_print(msg1, console = FALSE)
      log_quiet(msg1, msg = TRUE)
      message(msg1)
    } else {
      msg1 <- "Warning:"
      for (n in seq.int(to = 1, by = -1, length.out = sys.nframe() - 1)) {
        e1 <- sys.frame(n)
        # print(paste("frame: ", n))
        lse <- ls(e1)
        #print(lse)

        if ("call" %in% lse && "msg" %in% lse) {
          # call1 <- capture.output(print(get("call", envir = e1)))
          # print(msg)
          msg1 <- paste(msg1, get("msg", envir = e1))
          log_print(msg1, console = FALSE)
          logr:::log_quiet(msg1, msg = TRUE)
          message(msg1)
        }
      }
    }

    # Capture warnings locally
    # This is necessary because now the warnings() function in Base R
    # Doesn't work for logr.  So this allows a local version of the
    # warnings() function called get_warnings().
    wrn <- e$log_warnings
    wrn[length(wrn) + 1] <- paste0(msg1, collapse = "\n")
    e$log_warnings <- wrn

    # Publish warnings
    options("logr.warnings" = wrn)

  } else {

    logr:::disconnect_handlers()

  }
}
