#' @title Estimate Global Indices 
#' @description Estimate global indicators of the ecoystem state (richness, diversity, and other indicators).
#' @inheritParams global
#' @return A data.frame of samples x global indicators
#' @details This function returns the indices with the default choices for detection, prevalence and other parameters for simplicity and standardization. See the individual functions for more options. This function extends the functionality of the phyloseq::estimate_richness function.
#' @export
#' @seealso rarity, core_abundance, dominance, low_abundance, dominance, gini, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
diversity <- function(x) {

  .Deprecated(new = "global", msg = "The microbiome::diversity function has been replaced with the microbiome::global function and will be removed in the next release.")	  
  global(x)

}

