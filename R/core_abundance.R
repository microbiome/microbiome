#' @title Core Abundance Index
#' @description Calculates the community core abundance index.
#' @inheritParams core
#' @param split (Optional). Logical. Pool all samples and estimate index for the entire set.
#' @return A vector of core abundance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- core_abundance(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The core abundance index gives the relative proportion of the core species (in [0,1]). The core taxa are defined as those that exceed the given population prevalence threshold at the given detection level.
#' @seealso rarity
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection = .1/100, prevalence = 50/100, split = TRUE) {

  # Ensure the data is compositional	       
  xcomp <- transform(x, "compositional")	       

  # Core members
  cm <- core_members(xcomp, detection, prevalence)

  # Pick the core
  xc <- prune_taxa(cm, xcomp)

  colSums(abundances(xc))

}



