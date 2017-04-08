#' @title Estimate Global Indices 
#' @description Estimate global indicators of the ecoystem state (richness, diversity, and other indicators).
#' @param x \code{\link{phyloseq-class}} object
#' @param split (Optional). Logical. Should a separate set of richness
#'        estimates be performed for each sample? Or alternatively,
#'        pool all samples and estimate richness of the entire set.
#' @param measures Default is ‘NULL’, meaning that all available
#'          alpha-diversity measures will be included. Alternatively, you
#'          can specify one or more measures as a character vector of
#'          measure names. Values include those supported in the
#'          phyloseq::estimate_richness function:
#'          "Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher".
#'   In addition, the following measures are provided:
#'     "Richness" (number of unique taxa that give non-zero signal); 
#'     "Evenness" (Pielou's index);
#'     "Dominance" (Number of species needed to cover 50% of the ecosystem);
#'     "Top_Abundance" (Relative proportion of the most dominant species in [0,1]);
#'     "Rarity" (Relative proportion of the rare (non-core) species in [0,1]) - this complement (1-x) of the Core_Abundance
#'     "Low_Abundance" (Relative proportion of the least abundant species, below the detection level of 0.2%); 
#'     "Core_Abundance" (Relative proportion of the core species that exceed detection level 0.1% in over 50% of the samples);
#'     "Gini" (Gini index).
#' @inheritParams core
#' @return A data.frame of samples x global indicators; except when split=FALSE, a vector of indices is returned.
#' @details This function returns the indices with the default choices for detection, prevalence and other parameters for simplicity and standardization. See the individual functions for more options. This function extends the functionality of the phyloseq::estimate_richness function.
#' @examples 
#'   data(dietswap)
#'   d <- global(dietswap)
#' @export
#' @seealso rarity, core_abundance, top_abundance, low_abundance, dominance, gini, phyloseq::estimate_richness
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
diversity <- function(x, split = TRUE, measures = NULL) {

  .Deprecated(new = "global", msg = "The microbiome::diversity function has been replaced with the microbiome::global function and will be removed in the next release.")	  
  global(x, split = TRUE, measures = NULL)

}

