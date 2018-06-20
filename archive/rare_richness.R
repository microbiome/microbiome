#' @title Rare Richness
#' @description Calculate total rare richness
#' @param x \code{\link{phyloseq-class}} object
#' @param detection Detection threshold (non-negative real)
#' @param prevalence Prevalence threshold (in [0, 100])
#' @return Total rare abundance vector.
#' @examples
#'  data(dietswap)
#'  a <- rare_richness(dietswap)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rare_richness <- function(x, detection = 1, prevalence = 95) {
	      
  rare.taxa <- rare_members(x, detection = detection,
  	               prevalence = prevalence)

  # Rare matrix
  # FIXME: directly use the rare function
  xx <- abundances(prune_taxa(rare.taxa, x))
  xx <- matrix(xx, nrow = length(rare.taxa))

  # Rare taxa richness in each sample  
  ab <- colSums(xx > detection, na.rm = TRUE)

  ab

}
