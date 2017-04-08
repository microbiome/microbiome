#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @inheritParams core_abundance
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The rarity index gives the relative proportion of rare species (ie. those that are not part of the core microbiota) in the interval [0,1]. This is the complement (1-x) of the core abundance.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, rarity, diversity
rarity <- function(x, detection = .1/100, prevalence = 50/100, split = TRUE) {

  1 - core_abundance(x, detection, prevalence, split)

}

