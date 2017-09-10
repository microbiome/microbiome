#' @title Core Richness
#' @description Calculate total core richness
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold (non-negative real)
#' @param prevalence Prevalence threshold (in [0, 100])
#' @return Total core abundance vector.
#' @examples
#'  data(dietswap)
#'  a <- core_richness(dietswap)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_richness <- function(x, detection = 1, prevalence = 95) {
	      
  core.taxa <- core_members(x, detection = detection,
  	               prevalence = prevalence)

  # Core matrix
  # FIXME: directly use the core function
  xx <- abundances(prune_taxa(core.taxa, x))
  xx <- matrix(xx, nrow = length(core.taxa))

  # Core taxa richness in each sample  
  ab <- colSums(xx > detection, na.rm = TRUE)

  ab

}
