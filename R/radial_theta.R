#' @title Radial theta function.
#' @param x position parameter
#' @description Adapted from \pkg{NeatMap} and \pkg{phyloseq} packages but not exported and hence not available via phyloseq. Completely rewritten to avoid license conflicts. Vectorized to gain efficiency; only calculates theta and omits r. 
#' @keywords internal
radial_theta <- function (x) {

  x <- as(x, "matrix")
  theta <- atan2((x[, 2] - mean(x[, 2])), (x[, 1] - mean(x[, 1])))
  names(theta) <- rownames(x)
  
  theta
}