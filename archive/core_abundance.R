#' @title Core Abundance
#' @description Calculate total core abundance.
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold (non-negative real)
#' @param prevalence Prevalence threshold (in [0, 100])
#' @return Total core abunance vector.
#' @examples
#'  data(dietswap)
#'  a <- core_abundance(dietswap)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection = 1, prevalence = 95) {
	      
  core.taxa <- core_members(x, detection = detection,
  	               prevalence = prevalence)

  # Core matrix
  # FIXME: directly use the core function
  xx <- abundances(prune_taxa(core.taxa, x))
  xx <- matrix(xx, nrow = length(core.taxa))
  
  # Total sum of core abundances
  ab <- colSums(xx, na.rm = TRUE)

  ab

}
