#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @inheritParams core
#' @param split (Optional). Logical. Pool all samples and estimate index for the entire set.
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap, detection = 0.2/100)
#' @details The rarity index gives the relative proportion of rare species in [0,1]. The species that are below the indicated detection threshold are considered rare. Note that population prevalence is not considered.
#' @references
#'   Rarity has been used in this sense for instance in
#'   Kenneth J. Locey and Jay T. Lennon. Scaling laws predict global microbial diversity. PNAS 2016 113 (21) 5970-5975; doi:10.1073/pnas.1521291113.
#'
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
#' @seealso core_abundance, rarity, global
rarity <- function(x, detection = 0.2/100, split = TRUE) {

  # Ensure compositional data       
  xc <- transform(x, "compositional")

  if (!split) {
  
    otu <- rowSums(abundances(xc))
    x <- otu/sum(otu)
    do <- sum(x[x < detection])
    
  } else {

    otu <- abundances(xc)
    do <- apply(otu, 2, function (x) {sum(x[x < detection])})
    names(do) <- sample_names(x)

  }
  
  do

}



