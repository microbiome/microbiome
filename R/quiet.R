#' @title Quiet Output
#' @description Suppress all output from an expression. Works cross-platform.
#' @param expr Expression to run.
#' @param all If \code{TRUE} then suppress warnings and messages as well;
#' otherwise, only suppress printed output (such as from \code{print} or
#' \code{cat}).
#' @keywords internal
#' @return Used for its side effects.
#' @author
#' Adapted from \url{https://gist.github.com/daattali/6ab55aee6b50e8929d89}
#' @examples quiet(1 + 1)
#' @export
quiet <- function(expr, all=TRUE) {
    if (Sys.info()["sysname"] == "Windows") {
        file <- "NUL"
    } else {
        file <- "/dev/null"
    }
    
    if (all) {
        suppressWarnings(suppressMessages(
        suppressPackageStartupMessages(capture.output(expr, 
            file=file))))
    } else {
        capture.output(expr, file=file)
    }
}
