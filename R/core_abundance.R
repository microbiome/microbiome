#' @title Core Abundance
#' @description Calculate total core abundance
#' @param x phyloseq object
#' @param detection.threshold Detection threshold (non-negative real)
#' @param prevalence.threshold Prevalence threshold (in [0, 100])
#' @return Total core abunance vector.
#' @examples
#'  data(dietswap)
#'  a <- core_abundance(dietswap)
#' @export
#' @references See citation("microbiome")
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection.threshold = 1, prevalence.threshold = 95) {
	      
  core.taxa <- core(x, detection.threshold = detection.threshold,
  	               prevalence.threshold = prevalence.threshold)

  # Core matrix		       
  xx <- taxa_abundances(prune_taxa(core.taxa, x))

  # Total sum of core abundances
  ab <- colSums(xx, na.rm = TRUE)[rownames(xx)]

  ab

}
