#' @title Rarity Index
#' @description Calculates the community rarity index.
#' @inheritParams core
#' @param split (Optional). Logical. Should a separate set of richness
#'        estimates be performed for each sample? Or alternatively,
#'        pool all samples and estimate richness of the entire set.
#' @return A vector of rarity indices
#' @export
#' @examples
#'   data(dietswap)
#'   d <- rarity(dietswap, detection = 0.1/100, prevalence = 50/100)
#' @details The rarity index gives the relative proportion of rare species in [0,1]. 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
rarity <- function(x, detection, prevalence, split = TRUE) {

  # Ensure compositional data       
  xc <- transform(x, "compositional")

  # Pick rare members (those not in the core) 
  r <- taxa(rare(xc, detection = detection, prevalence = prevalence))

  if (!split) {
    otu <- abundances(xc)
    do <- sum((rowSums(otu)/sum(otu))[r])
  } else {

    otu <- abundances(xc)[r,]
    do <- colSums(otu)
    names(do) <- sample_names(x)
  }
  
  do

}



