#' @title Core Abundance Index
#' @description Calculates the community core abundance index.
#' @inheritParams rarity
#' @return A vector of core abundance indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- core_abundance(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The core abundance index gives the relative proportion of core species in [0,1]. This is exactly the complement (1 - x) of the rarity function; provided for convenience.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
core_abundance <- function(x, detection, prevalence, split = TRUE) {

  1 - rarity(x, detection, prevalence, split)

}



