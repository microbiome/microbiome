#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @inheritParams core_abundance
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The rarity index gives the relative proportion of rare species in [0,1]. The species that are below the indicated detection threshold are considered rare. Note that population prevalence is not considered.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, rarity, diversity
rarity <- function(x, detection = .1/100, prevalence = 50/100, split = TRUE) {

  1 - core_abundance(x, detection, prevalence, split)

}




