#' @title Core Taxa
#' @description Determine members of the core microbiota with given abundance and prevalence thresholds.
#' @param x phyloseq object
#' @param detection.threshold Detection threshold (non-negative real)
#' @param prevalence.threshold Prevalence threshold (in [0, 100])
#' @return Vector of core members
#' @examples
#'   data(dietswap)
#'   a <- core(dietswap, 1, 95)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core <- function(x, detection.threshold = 1, prevalence.threshold = 95)  {

  names(which(prevalence(x, detection.threshold) > prevalence.threshold))

}
