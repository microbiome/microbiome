#' @title Core Abundance Index
#' @description Calculates the community core abundance index.
#' @inheritParams core
#' @return A vector of core abundance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- core_abundance(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The core abundance index gives the relative proportion of the core species (in [0,1]). The core taxa are defined as those that exceed the given population prevalence threshold at the given detection level.
#' @seealso rarity
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection = .1/100, prevalence = 50/100) {

  # Ensure the data is compositional
  if (is.phyloseq(x)) {
    x <- abundances(x)
  }
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
  }
  # Compositional
  xcomp <- otu_table(apply(x, 2, function (x) {x/sum(x, na.rm = TRUE)})  , taxa_are_rows = TRUE)

  # Core members
  cm <- core_members(xcomp, detection, prevalence)

  # Pick the core
  xc <- prune_taxa(cm, xcomp)

  colSums(abundances(xc))

}



