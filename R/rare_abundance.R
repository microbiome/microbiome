#' @title Rare Abundance Index
#' @description Calculates the community rare_abundance index.
#' @inheritParams core_abundance
#' @return A vector of rare_abundance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rare_abundance(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The rare_abundance index gives the relative proportion of rare species (ie. those that are not part of the core microbiota) in the interval [0,1]. This is the complement (1-x) of the core abundance.
#' The rarity function provides the abundance of the least abundant taxa
#' within each sample, regardless of the population prevalence.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, rarity, diversity
rare_abundance <- function(x, detection = .1/100, prevalence = 50/100) {

  1 - core_abundance(x, detection, prevalence)

}

