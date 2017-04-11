#' @title Low Abundance Index
#' @description Calculates the concentration of low-abundance taxa below the indicated detection threshold.
#' @inheritParams core
#' @return A vector of indicators.
#' @export
#' @examples
#'   data(dietswap)
#'   d <- low_abundance(dietswap, detection = 0.2/100)
#' @details The low_abundance index gives the concentration of species at low abundance, or the  relative proportion of rare species in [0,1]. The species that are below the indicated detection threshold are considered rare. Note that population prevalence is not considered.
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, rarity, global
low_abundance <- function(x, detection = 0.2/100) {

  xc <- abundances(x, transform = "compositional")
  do <- apply(xc, 2, function (x) {sum(x[x < detection])})
  names(do) <- colnames(x)
  
  do

}



